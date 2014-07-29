/**@file  pgpr_chol.h
 * @brief This file provides Cholesky factorization and some related useful
 * functions such as inverse, log-determinant etc.
 *
 * @version 1.0
 */
#ifndef _PGPR_CHOL_H_
#define _PGPR_CHOL_H_
#include "pgpr_type.h"
/**@class pgpr_chol Cholesky factorization
 * @brief This class provides Cholesky factorization and some related useful functions
 * such as inverse, log-determinant etc.
 */
class pgpr_chol
{
private:
	Int n;   /**< The size of matrix, n x n */
	Mdoub el;  /**< The lower triangular Cholesky factor */
	Mdoub elt;  /**< The tranposed lower triangular matrix el */
public:
	/**@brief The constructor does a Cholesky factorization
	 *
	 * @detail The lower triangular factor is stored in el, and upper triangular factor
	 * is stored in elt.
	 *
	 * @param a the matrix that we want to do factorization
	 */
  pgpr_chol(Mdoub &A) :  el(A), elt(A) {
    Int i, j, k;
    Vdoub tmp;
    Doub sum;
		n = A.nrows();
    if (el.ncols() != n) throw("need square matrix");
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
        for (sum = el[i][j], k = i - 1; k >= 0; k--) sum -= el[i][k] * el[j][k];
        if (i == j) {
          if (sum <= 0.0)
            throw("pgpr_chol failed");
          el[i][i] = sqrt(sum);
        } else el[j][i] = sum / el[i][i];
      }
    }
    for (i = 0; i < n; i++) for (j = 0; j < i; j++) el[j][i] = 0.;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) elt[j][i] = el[i][j];
  }
	/**@brief Solving A* x = b using back substitution.
	 *
	 * @param b vector b = A * x
	 * @param x vector which would store the result.
	 */
  void solve(Vdoub &b, Vdoub &x) {
    Int i, k;
    Doub sum;
    if (b.size() != n || x.size() != n) throw("bad lengths in pgpr_chol");
    for (i = 0; i < n; i++) {
      for (sum = b[i], k = i - 1; k >= 0; k--) sum -= el[i][k] * x[k];
      x[i] = sum / el[i][i];
    }
    for (i = n - 1; i >= 0; i--) {
      for (sum = x[i], k = i + 1; k < n; k++) sum -= elt[i][k] * x[k];
      x[i] = sum / el[i][i];
    }
  }
  void elmult(Vdoub &y, Vdoub &b) {
    Int i, j;
    if (b.size() != n || y.size() != n) throw("bad lengths");
    for (i = 0; i < n; i++) {
      b[i] = 0.;
      for (j = 0; j <= i; j++) b[i] += el[i][j] * y[j];
    }
  }
  /** @brief solving el * y = b using back substitution.
   *
   *  @params b vector b = el * y
   *  @params y vector y stores the result
   */
  void elsolve(Vdoub &b, Vdoub &y) {
    Int i, j;
    Doub sum;
    if (b.size() != n || y.size() != n) throw("bad lengths");
    for (i = 0; i < n; i++) {
      for (sum = b[i], j = 0; j < i; j++) sum -= el[i][j] * y[j];
      y[i] = sum / el[i][i];
    }
  }

	///@brief Compute the inverse of a matrix
	///@param[out] ainv inverted matrix
  void inverse(Mdoub &ainv) {
    Int i, j, k;
    Doub sum;
    ainv.resize(n, n);
    for (i = 0; i < n; i++) for (j = 0; j <= i; j++) {
        sum = (i == j ? 1. : 0.);
        for (k = i - 1; k >= j; k--) sum -= el[i][k] * ainv[j][k];
        ainv[j][i] = sum / el[i][i];
      }
    for (i = n - 1; i >= 0; i--) for (j = 0; j <= i; j++) {
        sum = (i < j ? 0. : ainv[j][i]);
        for (k = i + 1; k < n; k++) sum -= elt[i][k] * ainv[j][k];
        ainv[i][j] = ainv[j][i] = sum / el[i][i];
      }
  }

  ///@brief Calculate the determinant of square matrix and retrun in log scale
	///@return Value of log-determinant
  Doub logdet() {
    Doub sum = 0.;
    for (Int i = 0; i < n; i++) sum += log(el[i][i]);
    return 2.*sum;
  }
};
#endif
