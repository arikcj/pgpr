/** @file pgpr_pic.h
 *  @brief This file provides the predictor (pgpr_pic)with PIC approximated Gaussian Process.
 */
#ifndef _PGPR_PIC_H_
#define _PGPR_PIC_H_
#include "pgpr_type.h"
#include "pgpr_cov.h"
#include "pgpr_chol.h"
#include "pgpr_cluster.h"
#include "pgpr_util.h"
//#define pic_tset_blk pic_tset_blk4
/** @class pgpr_pic
 *
 *  @brief This class provides the regression function using PIC Approximation.
 */
class pgpr_pic
{
private:
  pgpr_cov *cov;
  Doub h_mu;
  Vdoub pmu;
  Vdoub pvar;
  Doub elapsed;
  Doub rmse;
  Doub mnlp;

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

  Int pic_regr_blk( Mdoub data[], Vint ds,
                    Mdoub aset, Int as,
                    Mdoub tset[], Vint ts,
                    Vdoub pmu[], Vdoub pvar[]) {

    Int K = ds.size();
    Mdoub kuu;
    cov->se_ard(aset, as, kuu);
    pgpr_chol *chol_kuu = chol_cov(kuu);

    // local summary -> global summary & external summary
    Mdoub *ls_kuu;
    Vdoub *ls_zu;
    ls_kuu = new Mdoub[K];
    ls_zu = new Vdoub[K];
    Mdoub gs_kuu(kuu);
    Vdoub gs_zu(as, (Doub)0);
    Mdoub es_kuu(kuu);
    Vdoub es_zu(as, (Doub)0);

    pgpr_chol ** sdd_list = new pgpr_chol* [K];

    for (Int i = 0; i < K; i++) {
      //local summary in block i
      //pgpr_chol *chol_sdd = chol_pcov(data[i], ds[i], aset, as, chol_kuu);
      sdd_list[i] = chol_pcov(data[i], ds[i], aset, as, chol_kuu);
      pgpr_chol * chol_sdd = sdd_list[i];
      pitc_prep(data[i], ds[i], chol_sdd,
                aset, as, chol_kuu, ls_zu[i], ls_kuu[i]);
      //pitc_prep(data[i], ds[i], aset, as, chol_kuu,
      //          ls_zu, ls_kuu);
      //global summary
      for(Int r = 0; r < as; r++) {
        gs_zu[r] += ls_zu[i][r];
        for(Int c = 0; c < as; c++) {
          gs_kuu[r][c] += ls_kuu[i][r][c];
        }
      }
    }
    pgpr_chol *chol_suu = chol_cov(gs_kuu);
    //loop k
    for (Int k = 0; k < K; k++) {
      //external summary
      for(Int r = 0; r < as; r++) {
        es_zu[r] = gs_zu[r] - ls_zu[k][r];
        for(Int c = 0; c < as; c++) {
          es_kuu[r][c] = gs_kuu[r][c] - ls_kuu[k][r][c] - kuu[r][c];
        }
      }

      //predictive mean for each block
      Vdoub v(ds[k]);

      //pgpr_chol *chol_sdd = chol_pcov(data[k], ds[k], aset, as, chol_kuu);
      pgpr_chol *chol_sdd = sdd_list[k];

      Mdoub K_du;
      cov->se_ard(data[k], ds[k], aset, as, K_du);
      Mdoub K_td;
      cov->se_ard(tset[k], ts[k], data[k], ds[k], K_td);
      Mdoub K_tu;
      cov->se_ard(tset[k], ts[k], aset, as, K_tu);
#if 0
      Mdoub K_ut(as, ts[k]);

      for(int i = 0; i < ts[k]; i++) {
        for(int j = 0; j < as; j++) {
          K_ut[j][i] = K_tu[i][j];
        }
      }
#endif
      for(Int i = 0; i < ds[k]; i++) {
        v[i] = data[k][i][cov->dim] - h_mu;
      }
      if (pmu[k].size() != ts[k]) {
        pmu[k].resize(ts[k]);
      }
      Vdoub t1(ts[k]);
      Vdoub t2(ts[k]);
      Mdoub m_tu(ts[k], as);
      Mdoub tk_tu(ts[k], as);
      A_invB_C(K_td, chol_sdd, v, pmu[k]);
      A_invB_C(K_tu, chol_kuu, es_zu, t1);
      A_invB_C(K_tu, chol_kuu, es_kuu, m_tu);
      A_invB_C(K_td, chol_sdd, K_du, tk_tu);
      for(int i = 0; i < ts[k]; i++) {
        for(int j = 0; j < as; j++) {
          tk_tu[i][j] += m_tu[i][j];
        }
      }
      A_invB_C(tk_tu, chol_suu, gs_zu, t2);
      for(int i = 0; i < ts[k]; i++) {
        pmu[k][i] += (h_mu + t1[i] - t2[i]);
      }

      //predictive variance
      if (pvar[k].size() != ts[k]) {
        pvar[k].resize(ts[k]);
      }
      trace_A_invB_C(tk_tu, chol_suu, pvar[k]);
      trace_A_invB_transC(m_tu, chol_kuu, K_tu, t1);
      trace_A_invB_C(K_td, chol_sdd, t2);
      for(int i = 0; i < ts[k]; i++) {
        pvar[k][i] += cov->nos + cov->sig - t1[i] - t2[i];
      }
    }
    for( int i = 0; i < K; i++ ) {
      delete sdd_list[i];
    }
    delete [] sdd_list;
    return SUCC;
  }

  // pvar - the predictive variance
  Int pic_regr(  Mdoub data[], Vint ds,
                 Mdoub aset, Int as,
                 Mdoub tset, Int ts,
                 Vdoub &pmu, Vdoub &pvar) {
    pmsg(LEV_DBG, stdout, "PIC Regression\n");
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
    pgpr_cluster cluster(CLUSTER_ALGO2);
    cluster.pic_tset_blk(data, ds, tset, ts, tset_blk, ts_blk, pmu_blk, pvar_blk);

    pic_regr_blk(data, ds, aset, as, tset_blk, ts_blk, pmu_blk, pvar_blk);

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

    return SUCC;
  }

public:
  pgpr_pic(Char * hypf) {
    cov = new pgpr_cov(hypf);
    h_mu = cov->mu;
  }

  ~pgpr_pic() {
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
    pic_regr(blk_data, v_ds, supportset, ss, testset, ts,  pmu, pvar);
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
