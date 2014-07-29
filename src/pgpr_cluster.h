/** @file pgpr_cluster.h
 *  @brief This file provides some basic clustering methods for the GP model.
 */
#ifndef _PGPR_CLUSTER_H_
#define _PGPR_CLUSTER_H_
#include "pgpr_type.h"
#include "pgpr_cov.h"
#include "pgpr_chol.h"
#include "pgpr_util.h"
#define CLUSTER_ALGO0 0
#define CLUSTER_ALGO1 1
#define CLUSTER_ALGO2 2
#define CLUSTER_ALGO3 3
#define CLUSTER_ALGO4 4
/** @class pgpr_cluster
 *  @brief Provides some basic clustering algorithms used in the approximated GP algorithms.
 */
class pgpr_cluster
{
private:

  Int algo;

  void calculateMean( Mdoub data[], Vint ds, Mdoub* tset_blk, Vint & ts_blk, Mdoub & cents ) {
    int K = cents.nrows();
    int dimension = cents.ncols();
    for( int i = 0; i < K; i++ ) {
      for( int j = 1; j < ds[i]; j++ ) {
        for( int k = 0; k < dimension; k++ ) {
          cents[i][k] = data[i][j][k];
        }
      }
    }
    for( int i = 0; i < K; i++ ) {
      for( int j = 0; j < ts_blk[i]; j++ ) {
        for( int k = 0; k < dimension; k++ ) {
          cents[i][k] = tset_blk[i][j][k];
        }
      }
    }
    for( int i = 0; i < K; i++ ) {
      for( int k = 0; k < dimension; k++ ) {
        cents[i][k] /= (double) ( ds[i] + ts_blk[i] - 1 );
      }
    }
  }
  void distributeBlocks( Mdoub & cents, Mdoub& tset, Int ts, Mdoub* tset_blk, Vint & ts_blk ) {
    int K = cents.nrows();
    int dimension = cents.ncols();

    int ** allPointers = new int * [K];
    int *  pointSize =  new int[K];
    int *  inrSize =    new int[K];
    int blockSize = floor(ts / K);
    for( int i = 0; i < K; i++ ) {
      int bSize = blockSize;
      if( i == (K - 1) ) {
        bSize = ts - blockSize * ( K - 1 );
      }
      allPointers[i] = new int[ bSize ];
      pointSize[i] = bSize;
      inrSize [i] = 0;
    }

    for( int i = 0; i < ts; i++ ) {
      int minInd = -1;
      double minDis = DBL_MAX;
      for( int j = 0; j < K; j++ ) {
        double dis = 0;
        for( int k = 0; k < dimension; k++ ) {
          dis += (tset[i][k] - cents[j][k]) * (tset[i][k] - cents[j][k]);
        }
        if ( dis < minDis && inrSize[j] < pointSize[j]) {
          minInd = j ;
          minDis = dis;
        }
      }
      allPointers[minInd][inrSize[minInd]] = i;
      inrSize[minInd] ++ ;
    }

    for( int i = 0; i < K; i++ ) {
      for( int j = 0; j < pointSize[i]; j++ ) {
        int index = allPointers[i][j];
        takeSamp( tset[index], tset_blk[i][j], dimension );
        tset_blk[i][j][dimension] = index;
      }
    }

    for( int i = 0; i < K; i++ ) {
      delete [] allPointers[i];
    }
    delete [] allPointers;
    delete [] pointSize;
  }
  bool checkDistance( Mdoub & cents, Mdoub & cents2, double err ) {
    int K = cents.nrows();
    int dimension = cents2.ncols();
    double dis = 0;
    for( int i = 0; i < K; i++ ) {
      for( int j = 0; j < dimension; j++ ) {
        dis += (cents[i][j] - cents2[i][j] ) * (cents[i][j] - cents2[i][j] );
      }
    }
    dis = sqrt( dis );
    if( dis > err )
      return false;
    else
      return true;
  }


  /* Implement the clustering using the direct splitting
   */
  Int pic_tset_blk0( Mdoub data[], Vint ds, Mdoub tset, Int ts,
                     Mdoub* tset_blk, Vint & ts_blk, Vdoub* pmu_blk, Vdoub* pvar_blk ) {

    Int K = ds.size();
    Int bs = floor(ts / K);
    Int d_xz = tset.ncols();
    for (Int k = 0; k < K - 1; k++) {
      tset_blk[k].resize(bs, d_xz);
      ts_blk[k] = bs;
      pmu_blk[k].resize(bs);
      pvar_blk[k].resize(bs);

      for (Int i = 0; i < bs; i++) {
        Int pos = k * bs + i;
        takeSamp(tset[pos], tset_blk[k][i], d_xz);
      }
    }

    //The last block may contain slightly more points
    Int rs = ts - bs * (K - 1);
    tset_blk[K - 1].resize(rs, d_xz);
    ts_blk[K - 1] = rs;
    pmu_blk[K - 1].resize(rs);
    pvar_blk[K - 1].resize(rs);
    for (Int i = 0; i < rs; i++) {
      Int pos = bs * (K - 1) + i;
      takeSamp(tset[pos], tset_blk[K - 1][i], d_xz);
    }

    return SUCC;
  }
  /* Implement the clustering using random assign
   */
  Int pic_tset_blk1( Mdoub data[], Vint ds, Mdoub tset, Int ts,
                     Mdoub* tset_blk, Vint & ts_blk, Vdoub* pmu_blk, Vdoub* pvar_blk ) {

    Int K = ds.size();
    Int bs = floor(ts / K);
    Int d_xz = tset.ncols();
    int * list = new int[ts];
    for( int k = 0; k < ts; k++ ) {
      list[k] = k;
    }
    int numLeft = ts;
    for( Int k = 0; k < K - 1; k++ ) {
      tset_blk[k].resize( bs, d_xz );
      ts_blk[k] = bs;
      pmu_blk[k].resize( bs );
      pvar_blk[k].resize( bs );

      for( Int i = 0; i < bs; i++ ) {
        //Int pos = k * bs + i;
        int index = list[RANI(numLeft)];
        takeSamp( tset[index], tset_blk[k][i], d_xz );

        int t = list[numLeft - 1];
        list[index] = list[numLeft - 1];
        list[numLeft - 1] = t;
        numLeft --;
      }
    }
    Int rs = ts - bs * (K - 1);
    tset_blk[K - 1].resize(rs, d_xz);
    ts_blk[K - 1] = rs;
    pmu_blk[K - 1].resize(rs);
    pvar_blk[K - 1].resize(rs);
    for (Int i = 0; i < rs; i++) {
      takeSamp(tset[list[i]], tset_blk[K - 1][i], d_xz);
    }
    delete list;
    return SUCC;
  }

  Int pic_tset_blk2( Mdoub data[], Vint ds, Mdoub tset, Int ts,
                     Mdoub* tset_blk, Vint & ts_blk, Vdoub* pmu_blk, Vdoub* pvar_blk ) {
    //initalized data blocks
    Int K = ds.size();
    Int bs = ceil((Doub)ts / K);
    Int d_xz = tset.ncols();
    for (Int k = 0; k < K; k++) {
      tset_blk[k].resize(bs, d_xz);
      ts_blk[k] = 0;
    }
    //assign test points
    for (Int t = 0; t < ts; t++) {
      //shortest distance
      Doub sd = DBL_MAX;
      Int  ck = -1;
      Vdoub logdist(K);
      for (Int k = 0; k < K; k++) {
        Doub sq_dist = 0;
        for(Int i = 0; i < d_xz - 1; i++) {
          //the first point in data is the centroid
          sq_dist += SQR(data[k][0][i] - tset[t][i]);
        }
        sq_dist = SQRT(sq_dist);
        logdist[k] = sq_dist;
        if (sq_dist < sd && ts_blk[k] < bs) {
          ck = k;
          sd = sq_dist;
        }
      }
      //assign test point to closest block
      if(ck < 0) {
        throw("pic_tset_blk2: failed to assign point to a block\n");
      }
      Int cur = ts_blk[ck];
      takeSamp(tset[t], tset_blk[ck][cur], d_xz);
      //save the index of test set in measurement field
      tset_blk[ck][cur][d_xz - 1] = t;
      ts_blk[ck]++;
    }

    //resize pmu and pvar for each block
    for (Int k = 0; k < K; k++) {
      pmu_blk[k].resize(ts_blk[k]);
      pvar_blk[k].resize(ts_blk[k]);
    }
    return SUCC;
  }

  /* Implement the clustering using k-means
   */
  Int pic_tset_blk3( Mdoub data[], Vint ds, Mdoub tset, Int ts,
                     Mdoub* tset_blk, Vint & ts_blk, Vdoub* pmu_blk, Vdoub* pvar_blk ) {
    Int K = ds.size();
    Int blockSize = floor(ts / K);
    Int dimension = tset.ncols();
    dimension = dimension - 1; //don't count the measurement
    Mdoub cents( K, dimension);
    for( int i = 0; i < K; i++ ) {
      for( int j = 0; j < dimension; j++ ) {
        cents[i][j] = 0;
      }
    }

    for( int i = 0; i < K; i++ ) {
      for( int j = 0; j < ds[i]; j++ ) {
        for( int k = 0; k < dimension; k++ ) {
          cents[i][k] += data[i][j][k];
        }
      }
      for( int k = 0; k < dimension; k++ ) {
        cents[i][k] /= (double) ds[i];
      }
    }

    int ** allPointers = new int * [K];
    int *  pointSize =  new int[K];
    int *  inrSize =    new int[K];
    for( int i = 0; i < K; i++ ) {
      int bSize = blockSize;
      if( i == (K - 1) ) {
        bSize = ts - blockSize * ( K - 1 );
      }
      allPointers[i] = new int[ bSize ];
      pointSize[i] = bSize;
      inrSize [i] = 0;
    }

    for( int i = 0; i < ts; i++ ) {
      int minInd = -1;
      double minDis = DBL_MAX;
      for( int j = 0; j < K; j++ ) {
        double dis = 0;
        for( int k = 0; k < dimension; k++ ) {
          dis += (tset[i][k] - cents[j][k]) * (tset[i][k] - cents[j][k]);
        }
        if ( dis < minDis && inrSize[j] < pointSize[j]) {
          minInd = j ;
          minDis = dis;
        }
      }
      allPointers[minInd][inrSize[minInd]] = i;
      inrSize[minInd] ++ ;
    }

    for( int i = 0; i < K; i++ ) {
      tset_blk[i].resize( pointSize[i], dimension + 1 );
      ts_blk[i] = pointSize[i];
      pmu_blk[i].resize( pointSize[i] );
      pvar_blk[i].resize( pointSize[i] );
      for( int j = 0; j < pointSize[i]; j++ ) {
        int index = allPointers[i][j];
        takeSamp( tset[index], tset_blk[i][j], dimension + 1 );
        tset_blk[i][j][dimension] = index;
      }
    }

    for( int i = 0; i < K; i++ ) {
      delete [] allPointers[i];
    }
    delete [] allPointers;
    delete [] pointSize;
    return SUCC;
  }

  /* Implement the clustering using k-means
   */
  Int pic_tset_blk4( Mdoub data[], Vint ds, Mdoub tset, Int ts,
                     Mdoub* tset_blk, Vint & ts_blk, Vdoub* pmu_blk, Vdoub* pvar_blk ) {
    Int K = ds.size();
    Int blockSize = floor(ts / K);
    Int dimension = tset.ncols();
    dimension = dimension - 1; //don't count the measurement
    Mdoub cents( K, dimension );
    Mdoub cents2( K, dimension );
    for( int i = 0; i < K; i++ ) {
      for( int j = 0; j < dimension; j++ ) {
        cents[i][j] = 0;
        cents2[i][j] = 0;
      }
    }
    for( int i = 0; i < K; i++ ) {
      int bSize = blockSize;
      if( i == (K - 1) ) {
        bSize = ts - i * blockSize;
      }
      tset_blk[i].resize( bSize, dimension + 1 );
      ts_blk[i] = bSize;
    }

    for( int i = 0; i < K; i++ ) {
      for( int j = 0; j < ds[i]; j++ ) {
        for( int k = 0; k < dimension; k++ ) {
          cents[i][k] += data[i][j][k];
        }
      }
      for( int k = 0; k < dimension; k++ ) {
        cents[i][k] /= (double) ds[i];
      }
    }
    cents2 = cents;
    bool stop = false;
    const int ITER_MAX = 300000;
    int iter = 0;
    while( !stop && iter < ITER_MAX ) {
      distributeBlocks( cents, tset, ts, tset_blk, ts_blk );
      calculateMean( data, ds, tset_blk, ts_blk, cents2 );
      stop = checkDistance( cents, cents2, 0.000001 );
      if( stop == false ) {
        cents = cents2;
      }
      iter ++;
      //cout << "iter :" << iter << endl;
    }
    return SUCC;
  }

public:

  pgpr_cluster(Int al) {
    algo = al;
  }

  Int pic_tset_blk( Mdoub data[], Vint ds, Mdoub tset, Int ts,
                    Mdoub* tset_blk, Vint & ts_blk, Vdoub* pmu_blk, Vdoub* pvar_blk ) {
    switch(algo) {
    case CLUSTER_ALGO0:
      pic_tset_blk0( data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk );
      break;
    case CLUSTER_ALGO1:
      pic_tset_blk1( data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk );
      break;
    case CLUSTER_ALGO2:
      pic_tset_blk2( data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk );
      break;
    case CLUSTER_ALGO3:
      pic_tset_blk3( data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk );
      break;
    case CLUSTER_ALGO4:
      pic_tset_blk4( data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk );
      break;
    default:
      throw("The algorithm is not implemented\n");
      break;
    }
    return SUCC;
  }


  Int test_clustering(  Mdoub data[], Vint ds,
                        Mdoub aset, Int as,
                        Mdoub tset, Int ts,
                        Vdoub &pmu, Vdoub &pvar) {

    Int K = ds.size();
    Int bs = floor(ts / K);
    Int d_xz = tset.ncols();
    Mdoub *tset_blk;
    tset_blk = new Mdoub[K];
    Vint ts_blk(K);
    Vdoub *pmu_blk;
    pmu_blk = new Vdoub[K];
    Vdoub *pvar_blk;
    pvar_blk = new Vdoub[K];


    //Assign

    pic_tset_blk( data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk );

    return SUCC;
  }
};
#endif
