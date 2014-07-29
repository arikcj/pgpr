/** @file pgpr_pitc.h
 *  @brief This file provides the predictor with PITC (pgpr_pitc) Gaussian Process.
 */
#ifndef _PGPR_PITC_H_
#define _PGPR_PITC_H_
#include "pgpr_util.h"
#include "pgpr_cov.h"
#include "pgpr_chol.h"
/** @class pgpr_pitc
 *
 *  @brief This class provides the regression function using PITC Approximation.
 */
class pgpr_pitc
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

  Int pitc_regr_low3( Mdoub u, Int ss,
                      Mdoub suu, Vdoub fu,
                      Mdoub xt, Int ts, Int na,
                      Vdoub &t_mu, Vdoub &t_var) {

    Int K = na;
    Int bs = floor(ts / K);
    Int d_xz = xt.ncols();
    Mdoub *tset_blk;
    tset_blk = new Mdoub[K];
    Vint ts_blk(K);
    Vdoub *pmu_blk;
    pmu_blk = new Vdoub[K];
    Vdoub *pvar_blk;
    pvar_blk = new Vdoub[K];


    //Assign
    for (Int k = 0; k < K - 1; k++) {
      tset_blk[k].resize(bs, d_xz);
      ts_blk[k] = bs;
      pmu_blk[k].resize(bs);
      pvar_blk[k].resize(bs);

      for (Int i = 0; i < bs; i++) {
        Int pos = k * bs + i;
        takeSamp(xt[pos], tset_blk[k][i], d_xz);
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
      takeSamp(xt[pos], tset_blk[K - 1][i], d_xz);
    }

    //\Sigma_TT
    Vdoub alpha(ss);
    //Mdoub t_cov;
    Mdoub kuu;
    cov->se_ard(u, ss, kuu);
    //pmat(kuu);
    pgpr_chol chol_kuu(kuu);
    //the domain can be extremely large
    //Use post_var instead
    //  this matrix may require a large of memory
    pgpr_chol chol_suu(suu);

    //t_var.resize(ts);
    chol_suu.solve(fu, alpha);

    // tset_blk, ts_blk, pmu_blk, pvar_blk, alpha,
    // chol_suu, chol_kuu
    for( int i = 0; i < na; i++ ) {
      pitc_regr_low2_core( u, ss, tset_blk[i], ts_blk[i], pmu_blk[i], pvar_blk[i], alpha, chol_suu, chol_kuu);
    }

    //Use Mapreduce multicore


    /*
       for( int i = 0; i < na; i++ ) {
       post_var(u, ss, &chol_kuu, tset_blk[i], ts_blk[i], pvar_blk[i]);
       Mdoub K_td;
       cov->se_ard( tset_blk[i], ts_blk[i], u, ss, K_td);
       for( int k = 0; k < ts_blk[i]; k++ ) {
       pmu_blk[i][k] = h_mu;
       for(Int j = 0; j < ss; j++) {
       pmu_blk[i][k] += K_td[k][j] * alpha[j];
       }
       Vdoub K_ti(ss);
       for(Int j = 0; j < ss; j++) {
       K_ti[j] = K_td[k][j];
       }
       Vdoub v(ss);
       chol_suu.elsolve(K_ti, v);
       for(Int j = 0; j < ss; j++) {
       pvar_blk[i][k] += v[j] * v[j];
       }
       }
       }*/

    if (t_mu.size() != ts) {
      t_mu.resize(ts);
    }
    if( t_var.size() != ts ) {
      t_var.resize(ts);
    }
    for( int i = 0; i < na; i++ ) {
      for( int j = 0; j < ts_blk[i]; j++ ) {
        int ind = i * bs + j;
        t_mu[ind] = pmu_blk[i][j];
        t_var[ind] = pvar_blk[i][j];
      }
    }
    delete [] pmu_blk;
    delete [] pvar_blk;

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

  Int pitc_regr( Mdoub D[], Vint ds,
                 Mdoub aset, Int as,
                 Mdoub xt, Int ts,
                 Vdoub &t_mu, Vdoub &t_var) {


    Mdoub kuu;
    cov->se_ard(aset, as, kuu);
    pgpr_chol *chol_kuu = chol_cov(kuu);

    Mdoub ls_kuu;
    Vdoub ls_zu;
    Mdoub gs_kuu(kuu);
    Vdoub gs_zu(as, (Doub)0);
    for (Int i = 0; i < ds.size(); i++) {
      //local summary in block i
      pgpr_chol *chol_sdd = chol_pcov(D[i], ds[i], aset, as, chol_kuu);
      pitc_prep(D[i], ds[i], chol_sdd, aset, as, chol_kuu, ls_zu, ls_kuu);
      //global summary
      for(Int r = 0; r < as; r++) {
        gs_zu[r] += ls_zu[r];
        for(Int c = 0; c < as; c++) {
          gs_kuu[r][c] += ls_kuu[r][c];
        }
      }
    }
    pitc_regr_low3(aset, as, gs_kuu, gs_zu, xt, ts, ds.size(), t_mu, t_var);
    return SUCC;

  }

public:

  pgpr_pitc(Char * hypf) {
    cov = new pgpr_cov(hypf);
    h_mu = cov->mu;
  }

  ~pgpr_pitc() {
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
    pitc_regr(blk_data, v_ds, supportset, ss, testset, ts,  pmu, pvar);
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
