// Copyright (c) 2014, Jiangbo Yu

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file pgpr_plma.h
 *  @brief This file provides the main class of plma predictor
 *
 *  @author Jiangbo Yu
 */
#ifndef PGPR_PLMA_H_
#define PGPR_PLMA_H_
#include <Eigen/Dense>
#include "pgpr_parallel.h"
#include "mpi.h"
#include "pgpr_type.h"
#include "pgpr_chol.h"
#include "pgpr_cov.h"

#include <vector>
using namespace Eigen;
/**@class pgpr_plma  LMA parallel regression
 *  @brief This class provides the regression function using PLMA Approximation.
 */
class pgpr_plma
{
private:
	pgpr_parallel parallel;
	MatrixXd localTrainSamples;
	MatrixXd supportSamples;
	vector<VectorXd> *testClusters;
	pgpr_cov cov;
	int blocksize;
	int globalBlocks;
	int bandwidth;                // The number of samples should be read by current process
	int localband;
	double mean;                           // The mean
	double rmse;                           // The root mean square error
	double mnlp;                           
	MatrixXd pmu;                        // equal to a vector, predicted value
	int dim;                             // Dimension of the data
	MatrixXd pvar;                       // equal to a vector, predicted variance 
	MatrixXd trueval;                    // True value
	double elapsed;                        // Incurring time
	int myid;                            // The rank of the current process

	/**@brief Cluster the test data into training block,
	 *	   current implementation is greedy one.
	 *  @detail   The test points would belong to the nearest training block, where the distance is 
	 *     defined as the distance of  the current test point to the first train point in that
	 *     block. You should redefine a new clustering algorithm based on your own need. 
	 */
	
	void cluster( const MatrixXd & test){
		double sd = DBL_MAX;
		double sq_dist = 0;
		for(int testid = 0; testid < test.cols(); testid ++){
			// Compute the local minimal distance                                                                     
			//The first point in the data is the centroid of each cluster                                                    
			sd=0;
			for(int j = 0 ; j < dim; j ++){
				sd += SQR( localTrainSamples(j, 0) - test(j, testid));///cov->lsc[j]);      
			}
			sd = SQRT(sd);
			struct{
				double value;
				int virMachineID;
			} in, out;
			in.value = sd;
			in.virMachineID = myid;
			MPI_Barrier(MPI_COMM_WORLD);
			//use map-reduce method to locate the nearest cluster
			MPI::COMM_WORLD.Allreduce(&in, &out, 1, MPI::DOUBLE_INT, MPI::MINLOC);
			testClusters[out.virMachineID].push_back( test.col( testid));
		}
	}
public:
	/** @brief  Every machine loads the corresponding portion of data, and cluster the test data into nearest blocks.
	 *  @param hypf   the hyperparameter file name
	 *  @param train  file name of the training data
	 *  @param test   file name of the test data
	 *  @param supset file name of the support set
	 *  @param band   the bandwidth of LMA method
	 *  @param blks   the number of machines/blocks
	 */
     pgpr_plma(Char *hypf, Char *train, Char * test, Char * supset, int band, int blks)
		 :globalBlocks(blks), bandwidth(band), parallel(0), cov(hypf)
	{
		/*Be aware of that the Eigen is optimized using the column-major storeage
		 *The matrix should be read in a column-wise way.
		 */
		mean = cov.mu;
		dim = cov.dim;
		int dsize = getLines(train);
		int tsize = getLines(test);
		pmu = MatrixXd::Constant(tsize, 1, 0);
		pvar = MatrixXd::Constant(tsize, 1, 0);
		trueval = MatrixXd::Constant(tsize, 1, 0);
		parallel = pgpr_parallel(myid);
		blocksize = floor(dsize/globalBlocks);
		MPI_Comm_rank(MPI_COMM_WORLD,&myid);
		parallel.setRank(myid);
		//Load local train data
		localTrainSamples = MatrixXd(dim + 1, 1);
		/*When current machine is in the last serveral, bandwidth may not be as defined.
		 *For example, if band=1 && myid=2 && blk=3 then only one block in machine 3(myid=2)
		 *bandwidth should be redefined to be 0
		 */
		localband = MIN(bandwidth, globalBlocks - myid - 1);
		int ds = loadLocalData(train, localTrainSamples, blocksize, myid, localband + 1);
		
		//Load support sets
		supportSamples = MatrixXd(dim + 1, 1);
		int ss = loadData(supset, supportSamples);
		
		//Load test sets
		testClusters = new vector<VectorXd>[blks];
		MatrixXd testSamples(dim + 1, 1);
		loadData(test,testSamples);
		cluster(testSamples);
		int k = 0;
		for(int i = 0; i < globalBlocks; i ++)
			for(int j = 0; j < testClusters[i].size(); j ++)
				trueval(k++, 0) = testClusters[i][j][dim];
	}

	~pgpr_plma(){
		delete[] testClusters;	
	}
	
	/** @brief main funciton for LMA regression
	 *
	 * @detail The procedure is mainly following the AAAI2015 paper. Please read the paper
	 * to understand the full details of the implementation.
	 */
	void plma_regr(){

		//timer conters 
		pgpr_timer timer;
		
		//MPI
		MPI_Status status;
		MPI_Request send_request,recv_request;
		
		//K means kernel matrix, which is the same as Sigma used in the paper.
		// Compute \K__{S (D_m \cup D_m^B)} 
		timer.start();
		MatrixXd K_SS;
		Matrix<double, Dynamic, Dynamic, RowMajor> K_DS;
		MatrixXd K_DD;
		// K_DD_S = K_DD - K_DS * K_SS^-1 * K_SD
		MatrixXd K_DD_S;       
		MatrixXd tmp;
		//Check Noise
		cov.se_ard(supportSamples, K_SS);    
		cov.se_ard(localTrainSamples, supportSamples, tmp);
		K_DS = tmp;
		cov.se_ard_n(localTrainSamples, K_DD);

		MatrixXd K_DSTchol_SS; 
		MatrixXd K_SD;
		MatrixXd K_DmDm_SDmB;
		LLT<MatrixXd> lltofK_DmDm_SDmB;
		MatrixXd K_DmDmB_S;

		MatrixXd K_DmBDmB_S_inv;		
		MatrixXd dotK_SS;
		//MK_DmS = K_DmS - K_DmDmB_S * K_DmBDmB_S^-1 * K_DmB_S
		MatrixXd MK_DmS;

		
		/* --------------- Precompute redudant part ----------------------*/
		//        reS = K_DS * K_SS^-1
		MatrixXd reS;
		{
			MatrixXd K_SS_inv;
			K_SS_inv = K_SS.inverse();
			reS = K_DS * K_SS_inv;
		}
		/*-------------- Precompute end -------------------------------*/
		
		//TODO check noalias efficiency
		//K_DD_S = K_DD - K_DS * KSS^-1 * K_SD
		K_DD_S = K_DD;
		K_DD_S.noalias() -= reS * K_DS.transpose();
		
		MatrixXd reK;
		
	
        ///////////////////////////////
		//
		//   Compute dotK_SS
		//
		///////////////////////////////
		//Initialize: MK_DmS = K_DmS
		{
			/* --------------- Precompute redudant part ----------------------*/
			//        reK = K_DmDmB_S * K_DmBDmB_S^-1
			MatrixXd K_DmDmB_S;
			K_DmDm_SDmB = K_DD_S.topLeftCorner(blocksize, blocksize);
			if(localband >= 1){
				MatrixXd K_DmBDmB_S = K_DD_S.bottomRightCorner(localband * blocksize, localband * blocksize);
				K_DmDmB_S = K_DD_S.topRightCorner(blocksize, localband * blocksize);
				K_DmBDmB_S_inv = K_DmBDmB_S.inverse();
				reK = K_DmDmB_S * K_DmBDmB_S_inv;
			}
			/*-------------- Precompute end -------------------------------*/

			MK_DmS = K_DS.topRows(blocksize);
			
			//check localband when compute K_DmBDmB_SDmB
			if(localband >= 1){
				MatrixXd K_DmBS = K_DS.bottomRows(localband * blocksize);
				K_DmDm_SDmB.noalias() -= reK * K_DmDmB_S.transpose();
				//MK_DmS = K_DmS - K_DmDmB_S * K_DmBDmB_S^-1 * K_DmB_S
				MK_DmS.noalias() -= reK * K_DmBS;
			}
			
			// dotK_SS = MK_DmS^T * K_DmDm_SDmB^-1 * MK_DmS
			tmp = MK_DmS.transpose();
			lltofK_DmDm_SDmB = K_DmDm_SDmB.llt();
			lltofK_DmDm_SDmB.solveInPlace(MK_DmS);
			dotK_SS = tmp * MK_DmS;
			parallel.sync();
		}
		///////////////////////////////////
		//
		//    Compute dotY_S
		//
		////////////////////////////////////
			

		MatrixXd dotY_S;
		//MY_D = Y_D - K_DmDmB_S * KDmBDmB_S^-1 * YD
		VectorXd MY_Dm, Y_D;
		//Initialize MY_D = Y_D = y_Dm - mean
		Y_D = localTrainSamples.row(dim);
		Y_D.noalias() -= VectorXd::Constant(localTrainSamples.cols(), mean);
		MY_Dm = Y_D.head(blocksize);
		if(localband >= 1){
			VectorXd Y_DmB;
			Y_DmB = Y_D.tail(localband * blocksize);
			MY_Dm.noalias() -= reK * Y_DmB;
		}

		//dotY_S = MK_DmS^T * K_DmDm_SDmB^-1 * MY_D
		lltofK_DmDm_SDmB.solveInPlace(MY_Dm);
		dotY_S.noalias() = tmp * MY_Dm;
		
		//////////////////////////////////////
		//
		//     Compute barK_DU = barK_DU_S + Q_DU
		//     where Q_DU = K_DS * K_SS * K_SU
		//
		//////////////////////////////////////
		Matrix<double, Dynamic, Dynamic, RowMajor> *barK_DU = new Matrix<double, Dynamic, Dynamic, RowMajor>[globalBlocks];
		Matrix<double, Dynamic, Dynamic, RowMajor> *barK_DU_S = new Matrix<double, Dynamic, Dynamic, RowMajor>[globalBlocks];
		MatrixXd *Q_DU= new MatrixXd[globalBlocks];
		{
			for(int i = 0; i < globalBlocks; i ++)
				barK_DU_S[i].resize(localTrainSamples.cols(), testClusters[i].size()); 
			MatrixXd *K_SU = new MatrixXd[globalBlocks];
			for(int i = 0; i < globalBlocks; i ++){
				cov.se_ard(supportSamples, testClusters[i], tmp);
				K_SU[i] = tmp;
				Q_DU[i] = reS * K_SU[i];
			}
			// There exist some uncessary computing blocks here
			// ,which should not affact much. 
			for(int i = MAX(0, myid - bandwidth); i < MIN(myid + bandwidth + 1, globalBlocks); i ++){
				cov.se_ard(localTrainSamples, testClusters[i], tmp);
				barK_DU[i] = tmp;

				barK_DU_S[i] = barK_DU[i] - Q_DU[i];
			}
			
			delete[] K_SU;
		}

		/*-----------------------   Upper Triangular barK_DU_S------------------*/
		// It can be computed diagonally by diagonally.
		// Cur_diag is the global diagonal index working on.
		// We only need to start from bandwidth, the blocks inside the bandwidth
		// are already there. Note that the diag index starts from 0
		// when m + B < n: barK_DmUn_S = K_DmDmB_S * K_DmBDmB_S^-1 * barK_DmBUn_S
		int cur_diag = bandwidth + 1;
		while(cur_diag < globalBlocks){
		  
			/*--------------- Each machine computes local corresponding block---------*/
			
			// col is the global index of the test clusters/blocks. 
			int col = myid + cur_diag;
			
			if( col >= globalBlocks) break;
			if(bandwidth > 0){
				MatrixXd barK_DmBUn_S;
				barK_DmBUn_S = barK_DU_S[col].bottomRows(bandwidth * blocksize);
				tmp = reK * barK_DmBUn_S;
			}else{
				tmp = MatrixXd::Constant(blocksize, barK_DU_S[col].cols(), 0);
			}
			barK_DU_S[col].topRows(blocksize) = tmp;
			
			/*---------------- Communicate the sharing blocks -----------------*/
			//prepare for communication with the machine (myid - 1)
			int updateBlocks = MIN(cur_diag - bandwidth, bandwidth);
			if(myid > 0){
				tmp = barK_DU_S[col].topRows(updateBlocks * blocksize);
				parallel.send_msg(tmp, myid, myid - 1);
			}

			//receive from the machine myid + 1
			//some machine may be already idle now
			// if current machine havn't finish the last col, then the next machine would communicate with current one
			if(col < globalBlocks - 1 ){
				tmp.resize(updateBlocks * blocksize, testClusters[col + 1].size());
				parallel.recv_msg(tmp, myid + 1, myid + 1);
				barK_DU_S[col + 1].block(blocksize, 0, updateBlocks * blocksize, testClusters[col + 1].size()) = tmp;
			}
			cur_diag ++;
		}
		parallel.sync();

		/*-----------------------   Lower Triangular barK_DU_S------------------*/
		// Current implementation of computing lower triangular of barK_DU_S is by computing its 
		// tranpose barK_UD_S. In this way, we can save some computation and storage, but introduce
		// the communication overhead. It may be possible to be improved, depending on the communication
		// speed and memeory overhead. 
		// when m - B > n: barK_DmUn_S = (barK_UnDm_S)^T
		//                             = (barK_UnDnB_S * K_DnBDnB_S^-1 * barK_DnBDm_S)^T
		// Note barK_DnBDm_S is unknown, compute it using the similiar idea as computing upper triangular

		/* -------------- Precompute redudant part --------------*/
		//       reU = barK_UnDnB_S * K_DnBDnB_S^-1
		MatrixXd reU;
		if(localband >= 1){
			tmp = barK_DU_S[myid].bottomRows(localband * blocksize);
			reU = tmp.transpose() * K_DmBDmB_S_inv;
		}
		/*-------------- Precompute end --------------*/
		cur_diag = bandwidth + 1;
		
		Matrix<double, Dynamic, Dynamic, RowMajor> barK_DnBDm_S;
		if(localband == bandwidth)
			barK_DnBDm_S = K_DD_S.rightCols(blocksize);

		while( cur_diag < globalBlocks){
			int row, col;
			int updateBlocks = bandwidth;
			{   // Communication start: Sending the first bandwidth blocks of barK_DnDm_S to previous machine
				// @detail:  to compute barK_DmDn_S, machine m need barK_DmBDn_S which would be sent from machine m + 1
				col= myid + cur_diag - 1;
				if(col < globalBlocks){// other machines need do nothing
					if(myid > 0 && updateBlocks > 0){// sender 
						tmp = barK_DnBDm_S.topRows(updateBlocks * blocksize);
						parallel.send_msg(tmp, myid, myid - 1);
					}
					if(col  <  globalBlocks - 1 && updateBlocks > 0){// receiver
						tmp.resize(updateBlocks * blocksize, blocksize);
						parallel.recv_msg(tmp, myid + 1, myid + 1);
						barK_DnBDm_S.bottomRows(updateBlocks * blocksize) = tmp;
					}
				}
			}// Communication end;



			// Computer barK_DmUn_S^T  in machine n, which is barK_UnDm_S
			int upper_col = myid + cur_diag;
			if(upper_col < globalBlocks){
				if(bandwidth > 0){
					tmp = barK_DnBDm_S.bottomRows(bandwidth * blocksize);
					tmp = reU * tmp;
				}else{
					tmp = MatrixXd::Constant(testClusters[myid].size(), blocksize, 0);
				}

				// then send the result to the coressponding machine,which is transpose
				// sending the resulting barK_DmUn_S from machine n to machine m
				parallel.send_msg(tmp, myid, upper_col);
				if(bandwidth > 0)
					barK_DnBDm_S.block(0, 0, blocksize, blocksize) = reK * barK_DnBDm_S.bottomRows(bandwidth * blocksize);
			}

			//coressponding lower triangular col
			int lower_col = myid - cur_diag;
			if(lower_col>= 0){
				//the machines belong to the lower triangluar part should recive the message 
				tmp.resize(testClusters[lower_col].size(), blocksize);
				parallel.recv_msg(tmp, lower_col, lower_col);
				barK_DU_S[lower_col].topRows(blocksize) = tmp.transpose();
			}

			cur_diag ++;
			parallel.sync();
		}


		{//till now, each machine should have the complete first row of barK_DU_S
		
			// To fill the next localband rows of barK_DU_S, and reduce the communication overhead,
			// the messages are exchanged between only neighbouring machines.
			// for distance i, machine m send i + 1 row to machine m - 1.
			int distance = 1;
			while(distance <= bandwidth){

				double *sendbuf = NULL; 
				double *recvbuf = NULL;
				
				if(myid - bandwidth + distance > 1 && myid < globalBlocks - distance + 1){

					int cols = myid - bandwidth + distance - 1;
					int sendbufsize = 0;
					for(int i = 0; i < cols; i ++)
						sendbufsize += blocksize * testClusters[i].size();		     
					sendbuf = new double[sendbufsize];
					int k = 0;
					for(int col = 0; col < cols; col ++){
						for(int i = 0; i < blocksize; i ++)
							for(int j = 0; j < testClusters[col].size(); j ++)
								sendbuf[k++] = barK_DU_S[col](i + blocksize * (distance - 1), j);
					}
					MPI_Isend(sendbuf,sendbufsize, MPI_DOUBLE, myid - 1,myid, MPI_COMM_WORLD,&send_request);
					MPI_Wait(&send_request, &status);
				}
				
				if(myid + distance - bandwidth > 0 && myid < globalBlocks - distance ){
					//receiver
					int fromMachine = myid + 1;
					int recvbufsize = 0;
					int cols = myid + distance - bandwidth;
					for(int i = 0; i < cols; i ++)
						recvbufsize += blocksize * testClusters[i].size();
					recvbuf = new double[recvbufsize];
					
					MPI_Recv(recvbuf, recvbufsize, MPI_DOUBLE, fromMachine, fromMachine, MPI_COMM_WORLD, &status);
					
					int k = 0;
					for(int col = 0; col < cols ; col ++){
						for(int i = 0; i < blocksize; i ++)
							for(int j = 0; j < testClusters[col].size(); j ++)
								barK_DU_S[col](i + blocksize * distance, j) = recvbuf[k++];
					}
				}
				if(sendbuf != NULL){ delete sendbuf;sendbuf=0;}
				if(recvbuf != NULL) {delete recvbuf;recvbuf=0;}
				distance ++;
			}

		}
		parallel.sync();

		//    computer barK_DU = barK_DU_S + Q_DU
		for(int i = 0; i < globalBlocks; i ++)
			barK_DU[i] = barK_DU_S[i] + Q_DU[i];
		
		delete[] barK_DU_S;
		delete[]  Q_DU;

		VectorXd testOffset(globalBlocks);
		testOffset[0] = 0;
		for(int i = 1; i < globalBlocks; i ++)
			testOffset[i] = testOffset[i - 1] + testClusters[i - 1].size();


		//  Compute global summary ddotK_SS
		MatrixXd ddotK_SS;
		parallel.MapReduce(dotK_SS, ddotK_SS, MPI_SUM, -1);
		ddotK_SS.noalias() += K_SS;
		
		// Compute global summary ddotY_S
		MatrixXd ddotY_S;
		parallel.MapReduce(dotY_S, ddotY_S, MPI_SUM, -1);
		parallel.sync();
		
		//global summary
		MatrixXd ddotY_U;
		MatrixXd ddotK_US;
		MatrixXd ddotK_UU;
		// do prediction for each cluster/block
		/*  Here the prediciton is done in parallel. So each machine would store the coressponding part of dotY_U,
		 *  dotK_US and dotK_UU. For example, 
		 */
		for(int clusterid = 0; clusterid < globalBlocks; clusterid ++){
			int localtsize = testClusters[clusterid].size();
			if(localtsize <= 0 )
				continue;
			MatrixXd dotY_U = MatrixXd::Constant(localtsize, 1, 0);
			MatrixXd dotK_US = MatrixXd::Constant(localtsize, dotK_SS.rows(), 0);
			//Note: here we don't compute the whole matrix dotK_UU, but only diagonal part of dotK_UU.
			MatrixXd dotK_UU = VectorXd::Constant(localtsize, 1, 0);
			
			MatrixXd barK_DmU;
			barK_DmU = barK_DU[clusterid].topRows(blocksize);
			
			if(localband >= 1){
				barK_DmU.noalias() -= reK * barK_DU[clusterid].bottomRows(localband * blocksize);
			}
			tmp = barK_DmU.transpose();
			
			dotY_U.noalias() = tmp * MY_Dm;

			dotK_US.noalias() = tmp * MK_DmS;

			tmp = barK_DmU;
			lltofK_DmDm_SDmB.solveInPlace(barK_DmU);
			
			for(int i = 0; i < localtsize; i ++){
				dotK_UU(i, 0) = barK_DmU.col(i).dot(tmp.col(i));

			}
			parallel.MapReduce(dotY_U, ddotY_U, MPI_SUM, clusterid);
			parallel.MapReduce(dotK_US, ddotK_US, MPI_SUM, clusterid);
			parallel.MapReduce(dotK_UU, ddotK_UU, MPI_SUM, clusterid);
		}
		parallel.sync();
		delete[] barK_DU;
		//compute predicted mean and variance
		{
			//mean
			
			int offset = testOffset[myid];
			int localsize = testClusters[myid].size();
			MatrixXd localpmu = pmu;
			MatrixXd localpvar = pvar;
			if( localsize > 0){
				LLT<MatrixXd> lltofddotK_SS = ddotK_SS.llt();
				lltofddotK_SS.solveInPlace(ddotY_S);

				localpmu.block(offset, 0, localsize, 1).noalias() = ddotY_U + MatrixXd::Constant(localsize, 1, mean);
				localpmu.block(offset, 0, localsize, 1).noalias() -= ddotK_US * ddotY_S;
				
				//variance

				localpvar.block(offset, 0, localsize, 1)= MatrixXd::Constant(localsize, 1, cov.sig + cov.nos);
				tmp = ddotK_US.transpose();
				lltofddotK_SS.solveInPlace(tmp);
				
				for(int i = 0; i < localsize; i ++){
					localpvar(i + offset, 0) += ddotK_US.row(i).dot(tmp.col(i)) - ddotK_UU(i, 0);
				}
				
			}
			parallel.MapReduce(localpmu, pmu, MPI_SUM, 0);
			parallel.MapReduce(localpvar, pvar, MPI_SUM, 0);
		}
	}
	int regress(){
		pgpr_timer timer;
		parallel.sync(); timer.start();
		plma_regr();
		parallel.sync();
		elapsed = timer.end();
		if(myid == 0){
			rmse = getRmse(trueval.col(0), pmu.col(0));
			mnlp = getMnlp(trueval.col(0), pmu.col(0), pvar.col(0));
		}
		return 0;
	}

	void outputRst(Char * output){
		FILE *fp;
		int ts = pmu.size();
		fp = fopen(output, "w");
		if(fp == NULL){
			throw("Fail to open file\n");
		}
		for(int i = 0; i < ts; i++) {
			fprintf(fp,"%.4f %.4f %.4f\n",trueval(i,0), pmu(i, 0), pvar(i, 0));
		}
		pmsg(LEV_PRG, stdout, " %.4f | %.4f | %.4f |\n",elapsed,rmse, mnlp);

		fclose(fp);
	}
};
#endif
