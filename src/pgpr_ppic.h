/** @file pgpr_ppic.h
 *  @brief This file provides the predictor (pgpr_ppic) with parallel PIC Gaussian Process.
 */
#ifndef _PGPR_PPITC_H_
#define _PGPR_PPITC_H_
#include "mpi.h"
#include "pgpr_util.h"
#include "pgpr_cov.h"
#include "pgpr_chol.h"
#include "pgpr_cluster.h"
/** @struct t_local_summary
 *
 *  @brief This structure is used for local_summary.
 */
struct t_local_summary {
  Mdoub *ls_kuu;
  Vdoub *ls_zu;
} ;
/** @struct t_global_summary
 *
 *  @brief This structure is used for global_summary.
 */
struct t_global_summary {
  Mdoub *gs_kuu;
  Vdoub *gs_zu;
} ;
/** @class pgpr_ppic
 *
 *  @brief This class provides the regression function using PIC Approximation, but implemented in a paralle manner.
 */
class pgpr_ppic
{
private:
  pgpr_cov *cov;
  Doub h_mu;
  Vdoub pmu;
  Vdoub pvar;
  Doub elapsed;
  Doub rmse;
  Doub mnlp;

    static void addVec( Doub * ls, Doub * gs, Int * len, MPI_Datatype* t) {
    for( Int i = 0; i < *len; i++ ) {
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


  Int pic_core(Mdoub& datak, Int dsk,
               Mdoub & aset, Int as,
               Mdoub & tsetk, Int tsk,
               Vdoub & pmuk, Vdoub & pvark,
               Mdoub & kuu, pgpr_chol *chol_sddk,
               pgpr_chol * chol_kuu, pgpr_chol * chol_suu,
               t_local_summary & lvtk, t_global_summary & gvt) {
    Mdoub es_kuu(kuu);
    Vdoub es_zu(as, (Doub)0);
    for(Int r = 0; r < as; r++) {
      es_zu[r] = (*gvt.gs_zu)[r] - (*lvtk.ls_zu)[r];
      for(Int c = 0; c < as; c++) {
        es_kuu[r][c] = (*gvt.gs_kuu)[r][c] - (*lvtk.ls_kuu)[r][c] - kuu[r][c];
      }
    }

    //predictive mean for each block
    Vdoub v(dsk);

    //pgpr_chol *chol_sdd = chol_pcov(datak, dsk, aset, as, chol_kuu);

    Mdoub K_du;
    cov->se_ard(datak, dsk, aset, as, K_du);
    Mdoub K_td;
    cov->se_ard(tsetk, tsk, datak, dsk, K_td);
    Mdoub K_tu;
    cov->se_ard(tsetk, tsk, aset, as, K_tu);
    Mdoub K_ut(as, tsk);

    for(int i = 0; i < tsk; i++) {
      for(int j = 0; j < as; j++) {
        K_ut[j][i] = K_tu[i][j];
      }
    }

    for(Int i = 0; i < dsk; i++) {
      v[i] = datak[i][cov->dim] - h_mu;
    }
    if (pmuk.size() != tsk) {
      pmuk.resize(tsk);
    }

    Vdoub t1(tsk);
    Vdoub t2(tsk);
    Mdoub m_tu(tsk, as);
    Mdoub tk_tu(tsk, as);
    A_invB_C(K_td, chol_sddk, v, pmuk);
    A_invB_C(K_tu, chol_kuu, es_zu, t1);
    A_invB_C(K_tu, chol_kuu, es_kuu, m_tu);
    A_invB_C(K_td, chol_sddk, K_du, tk_tu);
    for(int i = 0; i < tsk; i++) {
      for(int j = 0; j < as; j++) {
        tk_tu[i][j] += m_tu[i][j];
      }
    }
    A_invB_C(tk_tu, chol_suu, (*gvt.gs_zu), t2);
    for(int i = 0; i < tsk; i++) {
      pmuk[i] += (h_mu + t1[i] - t2[i]);
    }

    //predictive variance
    if (pvark.size() != tsk) {
      pvark.resize(tsk);
    }
    trace_A_invB_C(tk_tu, chol_suu, pvark);
    trace_A_invB_C(m_tu, chol_kuu, K_ut, t1);
    trace_A_invB_C(K_td, chol_sddk, t2);
    for(int i = 0; i < tsk; i++) {
      pvark[i] += cov->nos + cov->sig - t1[i] - t2[i];
    }
    return SUCC;
  }

  Int pic_regr_blk_mpi( Mdoub data[], Vint ds,
                        Mdoub aset, Int as,
                        Mdoub tset[], Vint ts,
                        Vdoub pmu[], Vdoub pvar[]) {
    Mdoub kuu;
    cov->se_ard(aset, as, kuu);
    pgpr_chol *chol_kuu = chol_cov(kuu);

    Int K = ds.size();
    // local summary -> global summary & external summary
    int r;
    MPI_Comm_rank( MPI_COMM_WORLD, &r);

    Mdoub ls_kuu;
    Vdoub ls_zu;

    double * lu = new double [ as ];
    double * lk = new double [ as * as ];
    double * gu = new double [ as ];
    double * gk = new double [ as * as ];

    MPI_Op op1;
    MPI_Op_create( (MPI_User_function*)addVec, 1, &op1);

    Mdoub gs_kuu(kuu);
    Vdoub gs_zu(as, (Doub)0);
    pgpr_chol *chol_sdd = chol_pcov(data[r], ds[r], aset, as, chol_kuu);
    pitc_prep(data[r], ds[r], chol_sdd, aset, as, chol_kuu, ls_zu, ls_kuu);
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

    t_local_summary lvtk;
    lvtk.ls_kuu = &ls_kuu;
    lvtk.ls_zu = &ls_zu;
    t_global_summary gvt;
    gvt.gs_kuu = &gs_kuu;
    gvt.gs_zu = &gs_zu;
    pgpr_chol *chol_suu = chol_cov( gs_kuu );
    //loop k
    //

    pic_core( data[r], ds[r], aset, as, tset[r], ts[r], pmu[r], pvar[r], kuu, chol_sdd, chol_kuu, chol_suu, lvtk, gvt);

    int tsr = ts[r];

    double * pmur = new double [ tsr ];
    double * pvarr = new double [ tsr ];

    for( int i = 0; i < tsr; i++ ) {
      pmur[i] = pmu[r][i];
      pvarr[i] = pvar[r][i];
    }


    if( r == 0 ) {
      MPI_Status status;
      for( int i = 1; i < K; i++ ) {
        int tsi = ts[i];
        delete [] pmur;
        delete [] pvarr;
        pmur = new double[ tsi ];
        pvarr = new double[ tsi ];
        MPI_Recv(pmur, tsi, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status );
        MPI_Recv(pvarr, tsi, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status );

        pmu[i].resize( tsi );
        pvar[i].resize( tsi );

        for( int j = 0; j < tsi; j++ ) {
          pmu[i][j] = pmur[j];
          pvar[i][j] = pvarr[j];
        }
      }

    } else {
      MPI_Send ( pmur, tsr, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );
      MPI_Send ( pvarr, tsr, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD );
    }
    return SUCC;
  }

  Int ppic_regr (  Mdoub data[], Vint ds,
                   Mdoub aset, Int as,
                   Mdoub tset, Int ts,
                   Vdoub &pmu, Vdoub &pvar) {
#ifdef DEBUG
    pmsg(DBG, stdout, "PIC Regression\n");
#endif
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

    pgpr_cluster cluster(CLUSTER_ALGO2);
    cluster.pic_tset_blk( data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk );

    //pic regression
    pic_regr_blk_mpi(data, ds, aset, as, tset_blk, ts_blk, pmu_blk, pvar_blk);


    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank == 0 ) {
      //repack
      //Int cur = 0;
      for (Int k = 0; k < K; k++) {
        for (Int i = 0; i < ts_blk[k]; i++) {
          //pcp((Int)tset_blk[k][i][d_xz-1]);
          Int cur = (Int)tset_blk[k][i][d_xz - 1];
          pmu[cur] = pmu_blk[k][i];
          pvar[cur] = pvar_blk[k][i];
          //cur++;
        }
      }
    }

    return SUCC;
  }

public:
  pgpr_ppic(Char * hypf) {
    cov = new pgpr_cov(hypf);
    h_mu = cov->mu;
  }

  ~pgpr_ppic() {
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
    ppic_regr(blk_data, v_ds, supportset, ss, testset, ts,  pmu, pvar);
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
