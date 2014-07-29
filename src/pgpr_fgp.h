/** @file pgpr_fgp.h
 *  @brief This file provides the predictor (pgpr_fgp)with full Gaussian Process.
 */
#ifndef _PGPR_FGP_H_
#define _PGPR_FGP_H_
#include "pgpr_util.h"
#include "pgpr_cov.h"
#include "pgpr_chol.h"
/** @class pgpr_fgp
 *  @brief This class provides basic regression function using full GP.
 */
class pgpr_fgp
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

public:
  //full gp regression
  Int full_reg(Mdoub obs, Int dnum, Mdoub xt, Int ts, Vdoub &t_mu, Vdoub &t_var) {
    //step 1:cholesky decomposition
    Int ss = dnum;
    Vdoub v(ss);
    Vdoub alpha(ss);
    //Vdoub beta(ss);
    Mdoub K_dd;
    //build K+\sigma^2_n I.
    cov->se_ard_n(obs, ss, K_dd);
    pgpr_chol *chol = chol_cov(K_dd);

    //step 2: compute alpha
    for(Int i = 0; i < ss; i++) {
      v[i] = obs[i][cov->dim] - h_mu;
    }
    chol->solve(v, alpha);
    //step 3: predictive mean
    //Int ts = xt.nrows();
    t_mu.resize(ts);
    t_var.resize(ts);
    Mdoub K_td;
    cov->se_ard(xt, ts, obs, ss, K_td);

    // \bar{f}_*
    for(Int i = 0; i < ts; i++) {
      t_mu[i] = h_mu;
      for(Int j = 0; j < ss; j++) {
        t_mu[i] += K_td[i][j] * alpha[j];
      }

      //Step:4 predictive variance
      //k_{i,i}
      t_var[i] = cov->se_ard_n(xt[i], xt[i]);
      //k_{*f}
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
    }
    delete chol;
    return SUCC;
  }

  pgpr_fgp(Vdoub h, Int dx) {
    cov = new pgpr_cov(h, dx);
    h_mu = h[h.size() - 1];
  }

  pgpr_fgp(Char * hypf) {
    cov = new pgpr_cov(hypf);
    h_mu = cov->mu;
  }

  ~pgpr_fgp() {
    delete cov;
  }



  Int regress(Char * train, Char * test) {

    Int ddim = cov->dim + 1;
    Mdoub traindata(1, ddim);
    Mdoub testset(1, ddim);
    Int ds = loadData(train, traindata);
    Int ts = loadData(test, testset);
    pmu.resize(ts);
    pvar.resize(ts);

    pgpr_timer timer;
    timer.start();
    full_reg(traindata, ds, testset, ts, pmu, pvar);
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
