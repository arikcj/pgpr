/** @file pgpr_cov.h
 *  @brief This file provides the covariance class: pgpr_cov, which compute the
 *  covariance matrix.
 */
#ifndef _PGPR_COV_H_
#define _PGPR_COV_H_
#include "pgpr_util.h"
/** @class pgpr_cov
 *  @brief The pgpr_cov class provides the informaiton of covariance.
 */
 
class pgpr_cov
{
public:
	Doub nos;   /**< the noise variance of the data */	
	Vdoub lsc;  /**< vector of the length scale for each dimension */
	Doub sig;   /**< the signal variance of the data */
	Doub mu;    /**< the mean of the data */
	Int dim;    /**< the dimension of the features */
  /** @brief this functions initialized the class with a hyperparameter file
   *
   */
  pgpr_cov(Char * hypf) {
    FILE *fp = fopen(hypf, "r");
    if(fp == NULL) {
		throw("Fail to open file of hyperparameters\n");
    }
    fscanf(fp, "%lf ", &sig);
    fscanf(fp, "%lf ", &nos);
    fscanf(fp, "%lf ", &mu);
    fscanf(fp, "%d ", &dim);
    lsc.resize(dim);
    Doub tmp;
    for (Int i = 0; i < dim; i++) {
		fscanf(fp, "%lf ", &tmp);
		lsc[i] = tmp;
    }
    fprintf(fp, "\n");
    fclose(fp);
	sig = SQR(sig);
    nos = SQR(nos);

  }

  pgpr_cov(Vdoub h, Int d) {
    //h[]: noise,lenscale #1,...,lenscale #dim_x, signal, mean
    nos = SQR(h[0]);
    lsc.resize(d);
    for (Int i = 1; i <= d; i++) {
      lsc[i - 1] = h[i];
    }
    sig = SQR(h[d + 1]);
    dim = d;
  }

  inline Doub se_ard_n(const Doub *x, const Doub *y) {
    Doub val = 0;
    for (Int i = 0; i < dim; i++) {
      val += SQR((x[i] - y[i]) / lsc[i]);
    }

    return sig * exp(-0.5 * val) + nos;
  }

  inline Doub se_ard_n(const VectorXd &x, const VectorXd &y) {
    Doub val = 0;
    for (Int i = 0; i < dim; i++) {
      val += SQR((x[i] - y[i]) / lsc[i]);
    }

    return sig * exp(-0.5 * val) + nos;
  }
  inline Doub se_ard(const Doub *x, const Doub *y) {
    Doub val = 0;
    for (Int i = 0; i < dim; i++) {
      val += SQR((x[i] - y[i]) / lsc[i]);
    }

    return sig * exp(-0.5 * val);
  }

  inline Doub se_ard(const VectorXd &x, const VectorXd &y) {
    Doub val = 0;
    for (Int i = 0; i < dim; i++) {
      val += SQR((x[i] - y[i]) / lsc[i]);
    }

    return sig * exp(-0.5 * val);
  }
  inline void se_ard_n(Mdoub a, Int ss, Mdoub &k) {
    if (k.nrows() != ss) {
      k.resize(ss, ss);
    }
    for (Int i = 0; i < ss; i++) {
      for (Int j = i; j < ss; j++) {
        k[i][j] = se_ard(a[i], a[j]);
        if(i == j) {
          k[i][j] += nos;
        } else {
          k[j][i] = k[i][j];
        }
      }
    }
  }

  inline void se_ard_n(const MatrixXd &a, MatrixXd &k) {
	  Int ss = a.cols();
	  if (k.rows() != ss) {
		  k.resize(ss, ss);
	  }
	  for (Int i = 0; i < ss; i++) {
		  for (Int j = i; j < ss; j++) {
			  VectorXd r_i = a.col(i);
			  VectorXd r_j = a.col(j);
			  k(i,j) = se_ard(r_i, r_j);
			  if(i == j) {
				  k(i,j) += nos;
			  } else {
				  k(j,i) = k(i,j);
			  }
		  }
	  }
  }
  
  inline void se_ard(const Mdoub &a, Int ss, Mdoub &k) {
    if (k.nrows() != ss) {
      k.resize(ss, ss);
    }
    for (Int i = 0; i < ss; i++) {
      for (Int j = i; j < ss; j++) {
        k[i][j] = se_ard(a[i], a[j]);
        if(i != j) {
          k[j][i] = k[i][j];
        }
      }
    }
  }

  inline void se_ard(const MatrixXd &a, MatrixXd &k) {
	Int  ss = a.cols();
    if (k.rows() != ss) {
      k.resize(ss, ss);
    }
    for (Int i = 0; i < ss; i++) {
      for (Int j = i; j < ss; j++) {
		  VectorXd r_i = a.col(i);
		  VectorXd r_j = a.col(j);
		  k(i,j) = se_ard(r_i, r_j);
        if(i != j) {
			k(j,i) = k(i,j);
        }
      }
    }
  }
  

  //K_AB Mdoub version
  inline void se_ard(Mdoub a, Int ssa, Mdoub b, Int ssb, Mdoub &k) {
    k.resize(ssa, ssb);
    for (Int i = 0; i < ssa; i++) {
      for (Int j = 0; j < ssb; j++) {

        k[i][j] = se_ard(a[i], b[j]);
      }
    }

  }

  //K_AB Eigen version
  inline void se_ard(const MatrixXd &a, const MatrixXd &b, MatrixXd &k) {
	  Int ssa = a.cols();
	  Int ssb = b.cols();
	  k.resize(ssa, ssb);
	  for (Int i = 0; i < ssa; i++) {
		  for (Int j = 0; j < ssb; j++) {
			  VectorXd r_i = a.col(i);
			  VectorXd r_j = b.col(j);
			  k(i,j) = se_ard(r_i,r_j);
		  }
	  }
  }

//K_AB Eigen version
  inline void se_ard(const MatrixXd &a, const vector<VectorXd> &b, MatrixXd &k) {
	  Int ssa = a.cols();
	  Int ssb = b.size();
	  k.resize(ssa, ssb);
	  for (Int i = 0; i < ssa; i++) {
		  for (Int j = 0; j < ssb; j++) {
			  VectorXd r_i = a.col(i);
			  VectorXd r_j = b[j];
			  k(i,j) = se_ard(r_i,r_j);
		  }
	  }
  }
};
#endif
