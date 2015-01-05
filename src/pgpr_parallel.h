// Copyright (c) 2014, Jiangbo Yu

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/** @file pgpr_plma.h
 *  @brief This file provides the main class of plma predictor
 *
 *  @author Jiangbo Yu
 */
#ifndef PGPR_PARALLEL_H_
#define PGPR_PARALLEL_H_
#include <Eigen/Dense>
#include "mpi.h"
#include "pgpr_type.h"
#include "pgpr_chol.h"
#include "pgpr_cov.h"
#include <vector>
using namespace Eigen;
/**@class pgpr_parallel 
 *  @brief  the class provides the MPICH interface to commmunicate among the machines
 */
class pgpr_parallel
{
private:
	int rank;
public:
    pgpr_parallel(int r):rank(r){}
	void setRank(int r){ rank = r;}
/* @brief provide the interface to sending message using MPICH
	 * @param msg row-major matrix-structured message
	 * @param ori the message tag (the rank of orignation)
	 * @param dest the rank of destination 
	 */
	void send_msg(const Matrix<double, Dynamic, Dynamic, RowMajor> & msg, int ori, int dest){
		double *sendbuf = NULL;
		MPI_Request send_request;
		MPI_Status status;
		int bufsize;
		int rows = msg.rows();
		int cols = msg.cols();
		bufsize = rows * cols;
		sendbuf = new  double[bufsize];
		//row-major read row by row
		int k = 0, i = 0, j = 0;
		for(i = 0; i < rows ; i ++)
			for(j = 0; j < cols ; j ++)
				sendbuf[k++] = msg(i, j);
		MPI_Isend(sendbuf, bufsize, MPI_DOUBLE, dest, ori, MPI_COMM_WORLD, &send_request);
		MPI_Wait(&send_request, &status);
		delete sendbuf;
	}
	void send_msg(const MatrixXd &msg, int ori, int dest){
		double *sendbuf = NULL;
		MPI_Request send_request;
		MPI_Status status;
		int bufsize;
		int rows = msg.rows();
		int cols = msg.cols();
		bufsize = rows * cols;
		sendbuf = new  double[bufsize];
		//column-major read column by column
		int k = 0, i = 0, j = 0;
		for(j = 0; j < cols ; j ++)
			for(i = 0; i < rows ; i ++)
				sendbuf[k++] = msg(i, j);
		MPI_Isend(sendbuf, bufsize, MPI_DOUBLE, dest, ori, MPI_COMM_WORLD, &send_request);
		MPI_Wait(&send_request, &status);
		delete sendbuf;
	}

	/* @brief receive the message from the machine with rank ori
	 * @param msg the column-major matrix where the message would be stored
	 * @param ori the rank of origination
	 * @param flag the flag of the message
	 */
	void recv_msg(MatrixXd &msg, int ori, int flag){
		double *recvbuf = NULL;
		MPI_Status status;
		int bufsize;
		int rows = msg.rows();
		int cols = msg.cols();
		bufsize = rows * cols;
		recvbuf = new  double[bufsize];
		//column-major read column by column
		MPI_Recv(recvbuf, bufsize, MPI_DOUBLE, ori, flag, MPI_COMM_WORLD, &status);
		int k = 0, i = 0, j = 0;
		for(j = 0; j < cols ; j ++)
			for(i = 0; i < rows ; i ++)
				msg(i, j) = recvbuf[k++];
		delete recvbuf;
	}
	/* @brief provide the function to do global summary using mapich
	 * @param map the local matrix
	 * @param reduce the matrix used to store the result
	 * @param op MPICh summary operation 
	 * @param dest the machine id where reduced result should be stored, -1 means allreduce
	 
	 */
	void MapReduce(const MatrixXd &map, MatrixXd &reduce, MPI_Op op, int dest){
		double *sendbuf = NULL;
		double *recvbuf = NULL;
		int bufsize;
		int rows = map.rows();
		int cols = map.cols();
		bufsize = rows * cols;
		sendbuf = new  double[bufsize];
		recvbuf = new  double[bufsize];

		//column-major read column by column
		int k = 0, i = 0, j = 0;
		for(j = 0; j < cols ; j ++)
			for(i = 0; i < rows ; i ++)
				sendbuf[k++] = map(i, j);
		if(dest == -1)
			MPI_Allreduce(sendbuf, recvbuf, bufsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		else
			MPI_Reduce(sendbuf, recvbuf, bufsize, MPI_DOUBLE, op, dest, MPI_COMM_WORLD);
		if(dest == -1 || rank == dest){
			reduce.resize(rows, cols);
			k = 0;
			for(j = 0; j < cols ; j ++)
				for(i = 0; i < rows ; i ++)
					reduce(i, j) = recvbuf[k++];
		}
		delete sendbuf;
		delete recvbuf;
	}

	/* @brief provide the function to synchronize all the machinese. 
	 */
	void sync(){
		MPI_Barrier(MPI_COMM_WORLD);
	}
};
#endif
