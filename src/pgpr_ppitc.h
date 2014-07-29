/** @file pgpr_ppitc.h
 *  @brief This file provides the predictor (pgpr_ppitc) with parallel PITC Gaussian Process.
 */
#ifndef _PGPR_PPITC_H_
#define _PGPR_PPITC_H_
#include "mpi.h"
#include "pgpr_util.h"
#include "pgpr_cov.h"
#include "pgpr_chol.h"
/** @class pgpr_ppitc
 *
 *  @brief This class provides the regression function using PITC Approximation,implemented in a parallel manner.
 */
class pgpr_ppitc
{
private:
  pgpr_cov *cov;
  Doub h_mu;
  Vdoub pmu;
  Vdoub pvar;
  Doub elapsed;
  Doub rmse;
  Doub mnlp;

  static void addVec( double * ls, double * gs, int * len, MPI_Datatype* t) {
    for( int i = 0; i < *len; i++ ) {
      gs[i] += ls[i];
    }
  }

  inline pgpr_chol *chol_cov(Mdoub K_dd) {
    pgpr_chol *chol = new pgpr_chol(K_dd);
    return chol;
  }

  inline pgpr_chol *chol_cov(Mdoub obs, Int dnum) {
    Mdoub K_dd;
    cov->se_ard_n(obs, dnum, K_dd);
    pgpr_chol *chol = new pgpr_chol(K_dd);
    return chol;
  }

  inline pgpr_chol *chol_pcov(Mdoub obs, Int dnum, Mdoub act, Int anum) {
    Mdoub kdd;
    pgpr_chol *chol_kuu = chol_cov(act, anum);
    post_cov(act, anum, chol_kuu, obs, dnum, kdd);
    return chol_cov(kdd);
  }

  inline pgpr_chol *chol_pcov(Mdoub obs, Int dnum, Mdoub act, Int anum, pgpr_chol *chol_kuu) {
    Mdoub kdd;
    post_cov(act, anum, chol_kuu, obs, dnum, kdd);
    return chol_cov(kdd);
  }

  // @func  - compute posterior matrix using FGP
  Int post_var( Mdoub obs, Int ss, pgpr_chol *chol,
                Mdoub xt, Int ts,
                Vdoub &t_var) {

    //step 1:cholesky decomposition
    Vdoub v(ss);
    Vdoub beta(ss);
    //Mdoub K_dd;
    //build K+\sigma^2_n I.
    //cov->se_ard_n(obs,K_dd);
    //pgpr_chol chol(K_dd);

    //step 2: compute alpha

    //step 3: predictive mean
    if (t_var.size() != ts) {
      t_var.resize(ts);
    }
    Mdoub K_td;
    cov->se_ard(xt, ts, obs, ss, K_td);
    //pmat_s("ktd",K_td);

    //k_{T,T}
    //pmat_s("test",xt);
    //cov->se_ard_n(xt, ts, t_cov);


    // \bar{f}_*
    for(Int i = 0; i < ts; i++) {
      t_var[i] = cov->nos + cov->sig;

      //Step:4 predictive cov
      Vdoub K_ti(ss);
      for(Int j = 0; j < ss; j++) {
        K_ti[j] = K_td[i][j];
      }
      //var[f_ii]
      //v[i]=L\k*[i]
      chol->elsolve(K_ti, v);
      //predictive variance
      for(Int j = 0; j < ss; j++) {
        t_var[i] -= v[j] * v[j];
      }
      //var[\hat{f}_ji]
#if 0
      chol->solve(K_ti, beta);
      for(Int t = i + 1; t < ts; t++) {
        for(Int j = 0; j < ss; j++) {
          t_cov[t][i] -= K_td[t][j] * beta[j];
        }
        t_cov[i][t] = t_cov[t][i];
      }
#endif
    }
    return SUCC;

  }
  // @func  - compute posterior matrix using FGP
  Int post_cov( Mdoub obs, Int ss, pgpr_chol *chol,
                Mdoub xt, Int ts,
                Mdoub &t_cov) {
    //step 1:cholesky decomposition
    Vdoub v(ss);
    Vdoub beta(ss);
    //Mdoub K_dd;
    //build K+\sigma^2_n I.
    //cov->se_ard_n(obs,K_dd);
    //pgpr_chol chol(K_dd);

    //step 2: compute alpha

    //step 3: predictive mean
    t_cov.resize(ts, ts);
    Mdoub K_td;
    cov->se_ard(xt, ts, obs, ss, K_td);
    //pmat_s("ktd",K_td);

    //k_{T,T}
    //pmat_s("test",xt);
    cov->se_ard_n(xt, ts, t_cov);


    // \bar{f}_*
    for(Int i = 0; i < ts; i++) {

      //Step:4 predictive cov
      Vdoub K_ti(ss);
      for(Int j = 0; j < ss; j++) {
        K_ti[j] = K_td[i][j];
      }
      //var[f_ii]
      //v[i]=L\k*[i]
      chol->elsolve(K_ti, v);
      //predictive variance
      for(Int j = 0; j < ss; j++) {
        t_cov[i][i] -= v[j] * v[j];
      }
      //var[\hat{f}_ji]
      chol->solve(K_ti, beta);
      for(Int t = i + 1; t < ts; t++) {
        for(Int j = 0; j < ss; j++) {
          t_cov[t][i] -= K_td[t][j] * beta[j];
        }
        t_cov[i][t] = t_cov[t][i];
      }
    }
    return SUCC;
  }

  //prepare the local summary for PITC block
  Int pitc_prep( Mdoub D, Int ds, pgpr_chol *chol_sdd,
                 Mdoub U, Int us, pgpr_chol *chol_kuu,
                 Vdoub &fu, Mdoub &suu) {

    //\Sigma_{DD}|U+\sigma_n^2I
    //Mdoub sdd;
    //pcov(U, us, chol_kuu, D, ds, sdd);
    //pgpr_chol *chol_sdd = new pgpr_chol(sdd);
    Vdoub v(ds);
    Vdoub alpha(ds);
    Vdoub beta(ds);

    Mdoub K_ud;
    cov->se_ard(U, us, D, ds, K_ud);
    for(Int i = 0; i < ds; i++) {
      v[i] = D[i][cov->dim] - h_mu;
    }

    fu.resize(us);
    suu.resize(us, us);
    chol_sdd->solve(v, alpha);
    for(Int i = 0; i < us; i++) {
      //step 3: predictive mean
      fu[i] = 0;
      for(Int j = 0; j < ds; j++) {
        fu[i] += K_ud[i][j] * alpha[j];
      }
      //Step:4 predictive cov
      Vdoub K_ui(ds);
      for(Int j = 0; j < ds; j++) {
        K_ui[j] = K_ud[i][j];
      }
      //var[f_ii]
      //v[i]=L\k*[i]
      chol_sdd->elsolve(K_ui, v);
      suu[i][i] = 0;
      for(Int j = 0; j < ds; j++) {
        suu[i][i] += v[j] * v[j];
      }
      //var[\hat{f}_ij] where i!=j
      chol_sdd->solve(K_ui, beta);
      for(Int t = i + 1; t < us; t++) {
        suu[t][i] = 0;
        for(Int j = 0; j < ds; j++) {
          suu[t][i] += K_ud[t][j] * beta[j];
        }
        suu[i][t] = suu[t][i];
      }
    }
    return SUCC;
  }

  Int ppitc_regr(Mdoub D[], Vint ds,
                 Mdoub aset, Int as,
                 Mdoub xt, Int ts,
                 Vdoub &t_mu, Vdoub &t_var ) {
    Mdoub kuu;
    cov->se_ard(aset, as, kuu);
    pgpr_chol *chol_kuu = chol_cov(kuu);

    Mdoub ls_kuu;
    Vdoub ls_zu;
    Mdoub gs_kuu(kuu);
    Vdoub gs_zu(as, (Doub)0);
    MPI_Op op1;
    MPI_Op_create( (MPI_User_function*)addVec, 1, &op1);



    double * lu = new double[as];
    double * lk = new double[as * as];

    double * gu = new double[as];
    double * gk = new double[as * as];

    int r;
    MPI_Comm_rank( MPI_COMM_WORLD, &r);
    pgpr_chol *chol_sdd = chol_pcov(D[r], ds[r], aset, as, chol_kuu);
    pitc_prep(D[r], ds[r], chol_sdd, aset, as, chol_kuu, ls_zu, ls_kuu);


    for( int i = 0; i < as; i++ ) {
      lu[i] = ls_zu[i];
      gu[i] = 0;
      for( int j = 0; j < as; j++ ) {
        lk[i * as + j] = ls_kuu[i][j];
        gk[i * as + j] = 0;
      }
    }


    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( lu, gu, as, MPI_DOUBLE, op1, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Allreduce( lk, gk, as * as, MPI_DOUBLE, op1, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );
    for( int i = 0; i < as; i++ ) {
      gs_zu[i] = gu[i];
      for( int j = 0; j < as; j++ ) {
        gs_kuu[i][j] += gk[i * as + j];
      }
    }



    delete lu, delete [] lk, delete [] gu, delete [] gk;

    pitc_regr_low_mpi(aset, as, gs_kuu, gs_zu, xt, ts, ds.size(), t_mu, t_var);
    return SUCC;
  }


  Int pitc_regr_low_mpi( Mdoub u, Int ss,
                         Mdoub suu, Vdoub fu,
                         Mdoub xt, Int ts, Int na,
                         Vdoub &t_mu, Vdoub &t_var) {

    Int K = na;
    Int bs = floor(ts / K);
    Int d_xz = xt.ncols();

    Mdoub tset_blk;
    int ts_blk;
    Vdoub pmu_blk;
    Vdoub pvar_blk;

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if( rank < K - 1 ) {
      ts_blk = bs;
      tset_blk.resize( bs, d_xz );
      pmu_blk.resize( bs );
      pvar_blk.resize( bs );
      for (Int i = 0; i < bs; i++) {
        Int pos = rank * bs + i;
        takeSamp(xt[pos], tset_blk[i], d_xz);
      }
    } else {
      Int rs = ts - bs * (K - 1);
      ts_blk = rs;
      tset_blk.resize( rs, d_xz );
      pmu_blk.resize( rs );
      pvar_blk.resize( rs );
      for (Int i = 0; i < rs; i++) {
        Int pos = bs * (K - 1) + i;
        takeSamp(xt[pos], tset_blk[i], d_xz);
      }
      bs = rs;
    }
    //\Sigma_TT
    Vdoub alpha(ss);
    //Mdoub t_cov;
    Mdoub kuu;
    cov->se_ard(u, ss, kuu);
    //pmat(kuu)r
    pgpr_chol chol_kuu(kuu);
    //the domain can be extremely large
    //Use post_var instead
    //  this matrix may require a large of memory
    pgpr_chol chol_suu(suu);

    //t_var.resize(ts);
    chol_suu.solve(fu, alpha);

    // tset_blk, ts_blk, pmu_blk, pvar_blk, alpha,
    // chol_suu, chol_kuu
    pitc_regr_low2_core( u, ss, tset_blk, ts_blk, pmu_blk, pvar_blk, alpha, chol_suu, chol_kuu);

    //Use Mapreduce multicore

    if (t_mu.size() != ts) {
      t_mu.resize(ts);
    }
    if( t_var.size() != ts ) {
      t_var.resize(ts);
    }

    double * pmk = new double [bs];
    double * pvk = new double [bs];
    for( int i = 0; i < bs; i++ ) {
      pmk[i] = pmu_blk[i];
      pvk[i] = pvar_blk[i];
    }


    if( rank == 0 ) {
      MPI_Status status;
      for( int c = 0; c < bs; c++ ) {
        t_mu[c] = pmu_blk[c];
        t_var[c] = pvar_blk[c];
      }
      int recSize = 0;
      int start = 0;
      for( int i = 1; i < K; i++ ) {
        if( i < (K - 1) ) {
          recSize = ts / K;
        } else {
          recSize = ts - ts / K * i;
          delete [] pmk;
          delete [] pvk;
          pmk = new double[recSize];
          pvk = new double[recSize];
        }
        start = (ts / K) * i;
        MPI_Recv(pmk, recSize, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status );
        MPI_Recv(pvk, recSize, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status );
        for( int c = 0; c < recSize; c++ ) {
          t_mu[start + c] = pmk[c];
          t_var[start + c] = pvk[c];
        }
      }
    } else {
      MPI_Send(pmk, bs, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );
      MPI_Send(pvk, bs, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );
    }

    delete pmk, delete pvk;
    return SUCC;
  }

  void pitc_regr_low2_core( Mdoub& u, Int ss, Mdoub & tset_blki, int ts_blki, Vdoub& pmu_blki, Vdoub& pvar_blki, Vdoub& alpha, pgpr_chol& chol_suu, pgpr_chol & chol_kuu ) {
    post_var(u, ss, &chol_kuu, tset_blki, ts_blki, pvar_blki);
    Mdoub K_td;
    cov->se_ard( tset_blki, ts_blki, u, ss, K_td);
    for( int k = 0; k < ts_blki; k++ ) {
      pmu_blki[k] = h_mu;
      for(Int j = 0; j < ss; j++) {
        pmu_blki[k] += K_td[k][j] * alpha[j];
      }
      Vdoub K_ti(ss);
      for(Int j = 0; j < ss; j++) {
        K_ti[j] = K_td[k][j];
      }
      Vdoub v(ss);
      chol_suu.elsolve(K_ti, v);
      for(Int j = 0; j < ss; j++) {
        pvar_blki[k] += v[j] * v[j];
      }
    }
  }

public:
  pgpr_ppitc(Char * hypf) {
    cov = new pgpr_cov(hypf);
    h_mu = cov->mu;
  }
 
  ~pgpr_ppitc() {
    delete cov;
  }

  Int regress(Char * train, Char * test, Char * support, Int blocks) {
    Int ddim = cov->dim + 1;
    Mdoub traindata(1, ddim);
    Mdoub testset(1, ddim);
    Mdoub supportset(1, ddim);
    Int ds = loadData(train, traindata);
    //data
    Mdoub *blk_data = new Mdoub[blocks];
    Vint v_ds(blocks);
    Int pos = 0;
    Int blksize = floor(ds / blocks);
    Int remain = ds % blocks;

    // First Block
    blk_data[0].resize(blksize + remain, ddim);
    for (Int j = 0; j < blksize + remain; j++) {
      takeSamp(traindata[pos], blk_data[0][j], ddim);
      pos++;
    }
    v_ds[0] = blksize + remain;

    for(Int i = 1; i < blocks; i++) {
      blk_data[i].resize(blksize, ddim);
      for (Int j = 0; j < blksize; j++) {
        takeSamp(traindata[pos], blk_data[i][j], ddim);
        pos++;
      }
      v_ds[i] = blksize;
    }
    Int ts = loadData(test, testset);
    Int ss = loadData(support, supportset);
    pmu.resize(ts);
    pvar.resize(ts);
    pgpr_timer timer;
    timer.start();
    ppitc_regr(blk_data, v_ds, supportset, ss, testset, ts,  pmu, pvar);
    elapsed = timer.end();
    Vdoub trueval(ts);
    for(Int i = 0; i < ts; i++) {
      trueval[i] = testset[i][ddim - 1];
    }
    rmse = getRmse(trueval, pmu);
    mnlp = getMnlp(trueval, pmu, pvar);
    return SUCC;
  }
 
  void outputRst(Char * output) {
    FILE * fp;
    Int ts = pmu.size();
    fp = fopen(output, "w");
    if(fp == NULL) {
      throw("Fail to open file\n");
    }
    for(Int i = 0; i < ts; i++) {
      fprintf(fp, "%.4f %.4f\n", pmu[i], pvar[i]);
    }
    pmsg(LEV_PRG, stdout, " %.4f | %.4f | %.4f |\n", elapsed, rmse, mnlp);

    fclose(fp);

  }



};
#endif
