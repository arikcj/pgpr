/**
 * @file pgpr_util.h
 *	@brief This file contains a collection of useful functions such as
 *	macro-like inline functions, debug functions, exception handling,
 *	file I/O, and a real-time timer class. 
 *
 * @author CHEN jie, arik.cj@gmail.com
 * 
 * @version 1.0
 **/
#ifndef _PGPR_UTIL_H_
#define _PGPR_UTIL_H_
#include <sys/time.h>
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif
#include "pgpr_type.h"
#include "pgpr_chol.h"

#define LEV_PRG 0
#define LEV_DBG 5
#define LEV_ALL 10
#define CURLEV LEV_PRG

// macro-like inline functions

template<class T>
inline T SQR(const T a)
{
  return a * a;
}

template<class T>
inline const T &MAX(const T &a, const T &b)
{
  return b > a ? (b) : (a);
}

inline float MAX(const double &a, const float &b)
{
  return b > a ? (b) : float(a);
}

inline float MAX(const float &a, const double &b)
{
  return b > a ? float(b) : (a);
}

template<class T>
inline const T &MIN(const T &a, const T &b)
{
  return b < a ? (b) : (a);
}

inline float MIN(const double &a, const float &b)
{
  return b < a ? (b) : float(a);
}

inline float MIN(const float &a, const double &b)
{
  return b < a ? float(b) : (a);
}

template<class T>
inline T SIGN(const T &a, const T &b)
{
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const float &a, const double &b)
{
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float SIGN(const double &a, const float &b)
{
  return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

template<class T>
inline void SWAP(T &a, T &b)
{
  T dum = a;
  a = b;
  b = dum;
}

#define RANI(x) (Int)((Doub)rand()/RAND_MAX*x)
#define SRAND(x) srand(x); rand()
#define LOG(x) log(x)
#define SQRT(x) sqrt(x)
#define EXP(x) exp(x)
#define POW(x,y) (Int)pow((Doub)(x),(Doub)(y))

#define takeSamp(a,b,d) for(int _i=0; _i<d; _i++) b[_i]=a[_i]

// file I/O
inline Int getLines(Char * file)
{
  Int number_of_lines = 0;
  FILE *infile = fopen(file, "r");
  Int ch;

  while (EOF != (ch = getc(infile)))
    if ('\n' == ch)
      ++number_of_lines;
  return number_of_lines;
}

inline Int saveData(Char * file, Mdoub m_data)
{
  FILE * fp = fopen(file, "w");
  if(fp == NULL) {
    throw("save data: Fail to open file\n");
  }
  Int ddim = m_data.ncols();
  Int dsize = m_data.nrows();
  for(Int i = 0; i < dsize; i++) {
    fprintf(fp, "%.4f ", m_data[i][ddim - 1]);
    for(Int j = 0; j < ddim - 1; j++) {
      fprintf(fp, "%d:%.4f ", j, m_data[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return SUCC;
}

inline Int loadData(Char * file, Mdoub &m_data)
{
  FILE * fp = fopen(file, "r");
  if(fp == NULL) {
    throw("load data: Fail to open file\n");
  }
  Int ddim = m_data.ncols();
  Int dsize = getLines(file);
  m_data.resize(dsize, ddim);
  Doub tmp;
  Int i_tmp;
  for(Int i = 0; i < dsize; i++) {
    fscanf(fp, "%lf ", &tmp);
    m_data[i][ddim - 1] = tmp;
    for(Int j = 0; j < ddim - 1; j++) {
      fscanf(fp, "%d:%lf ", &i_tmp, &tmp);
      m_data[i][j] = tmp;
    }
    fscanf(fp, "\n");
  }
  fclose(fp);
  return dsize;
}

inline Int saveHyper(Char * file, Vdoub h, Int d)
{
  Doub nos = h[0];
  Vdoub lsc(d);
  for (Int i = 1; i <= d; i++) {
    lsc[i - 1] = h[i];
  }
  Doub sig = h[d + 1];
  Doub h_mu = h[h.size() - 1];
  //save Hyperparameter
  FILE *fp = fopen(file, "w");
  if(fp == NULL) {
    throw("Fail to open file for storing hyperparameters\n");
  }

  fprintf(fp, "%.4f ", sig);
  fprintf(fp, "%.4f ", nos);
  fprintf(fp, "%.4f ", h_mu);
  Int nlsc = lsc.size();
  fprintf(fp, "%d ", nlsc);
  for (Int i = 0; i < nlsc; i++) {
    fprintf(fp, "%.4f ", lsc[i]);
  }
  fprintf(fp, "\n");

  fclose(fp);
  return SUCC;
}

//Debug

#if defined(NODEBUG)
#define pmsg(level, outfp, format, args...) ((void)0)
#else
inline void pmsg(int level, FILE *outfp, const char* format, ...)
{
  va_list args;

  if (level > CURLEV)
    return;

  va_start(args, format);
  vfprintf(outfp, format, args);
  va_end(args);
}
/* print a message, if it is considered significant enough.
      Adapted from [K&R2], p. 174 */
#endif

#define pcp(x) pmsg(LEV_DBG,stdout,"CheckPoint> %d\n",x)
#define pcpst(x) pmsg(LEV_DBG,stdout,"[%d] start>\n",x)
#define pcpen(x) pmsg(LEV_DBG,stdout,"[%d] end# \n",x)
#define pcp_s(s,x) pmsg(LEV_DBG,stdout," %s > %.4lf\n",s,(double)x)

#define pvec(vec)\
    for(Int _j=0; _j<vec.size(); _j++)\
    {\
        pmsg(LEV_DBG,stdout,"%.4f ",(double)vec[_j]);\
    }\
    pmsg(LEV_DBG,stdout,"\n");

#define pvec_r(vec,vsize)\
    for(Int _j=0; _j<vsize; _j++)\
    {\
        pmsg(LEV_DBG,stdout,"%.4f ",(double)vec[_j]);\
    }\
    pmsg(LEV_DBG,stdout,"\n");


#define pvec_s(_s,vec)\
    pmsg(LEV_DBG,stdout,"%s:\n",_s);\
    for(Int _j=0; _j<vec.size(); _j++)\
    {\
        pmsg(LEV_DBG,stdout,"%.4f ",(double)vec[_j]);\
    }\
    pmsg(LEV_DBG,stdout,"\n");


#define pmat_r(mat,_r)\
    for(Int _i=0; _i<_r; _i++)\
    {\
        for(Int _j=0; _j<_r; _j++)\
        {\
            pmsg(LEV_DBG,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LEV_DBG,stdout,"\n");\
    }\
 
#define pmat_pos(mat)\
    for(Int _i=0; _i<mat.nrows(); _i++)\
    {\
        for(Int _j=0; _j<mat.ncols(); _j++)\
        {\
            if(mat[_i][_j]<0) break;\
            pmsg(LEV_DBG,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LEV_DBG,stdout,"\n");\
    }\
 

#define pmat(mat)\
    for(Int _i=0; _i<mat.nrows(); _i++)\
    {\
        for(Int _j=0; _j<mat.ncols(); _j++)\
        {\
            pmsg(LEV_DBG,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LEV_DBG,stdout,"\n");\
    }\
 
#define pmat_s(_s,mat)\
    pmsg(LEV_DBG,stdout,"%s:\n",_s);\
    for(Int _i=0; _i<mat.nrows(); _i++)\
    {\
        for(Int _j=0; _j<mat.ncols(); _j++)\
        {\
            pmsg(LEV_DBG,stdout,"%.4f ",(double)mat[_i][_j]);\
        }\
        pmsg(LEV_DBG,stdout,"\n");\
    }\
 
// Exception handling from Numerical Recipe
// usage example:
//
//	try {
//		somebadroutine();
//	}
//	catch(NRerror s) {NRcatch(s);}
//
// (You can of course substitute any other catch body for NRcatch(s).)

#ifndef _USENRERRORCLASS_
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#else
struct NRerror {
  char *message;
  char *file;
  int line;
  NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch(NRerror err)
{
  printf("ERROR: %s\n     in file %s at line %d\n",
         err.message, err.file, err.line);
  exit(1);
}
#endif


//Useful routines

//bubble sort
inline void bubble_sort(Vdoub a, Vint &ai) {
    Int j = 1;
    Int k = 1;
    Doub t = 1.0;
    Int n = a.size();
    while (n-- && t)
      for (t = j = k = 0; j < n; j++) {
        if (a[j] <= a[j + 1]) continue;
        t = a[j], a[j] = a[j + 1], a[j + 1] = t;
        t = ai[j], ai[j] = ai[j + 1], ai[j + 1] = t;
        t = 1;
      }
  }

// argmax_i v(i)
inline Int argmaxi(Vdoub v)
{
  Doub mval = 0;
  Int mi = -1;
  Int vs = v.size();
  Vint vm(vs);
  Int vcur = 0;
  for(Int i = 0; i < vs; i++) {
    if (mval < v[i]) {
      vcur = 0;
      mval = v[i];
      vm[vcur] = i;
      vcur++;
    } else if (mval == v[i]) {
      mval = v[i];
      vm[vcur] = i;
      vcur++;
    }
  }
  if (vcur > 0) {
    mi = vm[RANI(vcur)];
  }
  return mi;
}

// Root mean square error (see section 6.1)
inline Doub getRmse(Vdoub v1, Vdoub v2)
{
  Int sz = v1.size();
  Doub rmse = 0;
  if (sz == v2.size()) {
    for(Int i = 0; i < sz; i++) {
      rmse += SQR(v1[i] - v2[i]);
    }
    rmse = SQRT(rmse / sz);
  } else {
    throw("The two vectors have different size\n");
  }
  return rmse;
}

// Mean Negative Log Probability (see section 6.1)
inline Doub getMnlp(Vdoub v1, Vdoub v2, Vdoub v3)
{
  Int sz = v1.size();
  Doub mnlp = 0;
  if (sz == v2.size()) {
    for(Int i = 0; i < sz; i++) {
      mnlp += SQR(v1[i] - v2[i]) / v3[i] +  LOG(2 * 3.1416 * v3[i]);
    }
    mnlp = 0.5 * mnlp / sz;
  } else {
    throw("The two vectors have different size\n");
  }
  return mnlp;
}


//A: m by n matrix
//B: n by n matrix
//c: n-size vector
//D = A B^{-1} c
inline Int A_invB_C(Mdoub A, pgpr_chol *chol_b, Vdoub C, Vdoub &D)
{
  Int m = A.nrows();
  Int n = A.ncols();
  Vdoub beta(n);
  D.assign(m, 0);
  chol_b->solve(C, beta);
  for(Int k = 0; k < m; k++) {
    for(Int j = 0; j < n; j++) {
      D[k] += A[k][j] * beta[j];
    }
  }
  return SUCC;
}


//A: m by n matrix
//B: n by n matrix
//C: n by l matrix
//D = A B^{-1} C
inline Int A_invB_C(Mdoub A, pgpr_chol *chol_b, Mdoub C, Mdoub &D)
{
  Int m = A.nrows();
  Int n = A.ncols();
  Int l = C.ncols();
  Vdoub v(n);
  Vdoub beta(n);
  Vdoub Cn(n);
  D.assign(m, l, 0);
  for(Int i = 0; i < l; i++) {
    for(Int j = 0; j < n; j++) {
      Cn[j] = C[j][i];
    }
    chol_b->solve(Cn, beta);
    for(Int k = 0; k < m; k++) {
      for(Int j = 0; j < n; j++) {
        D[k][i] += A[k][j] * beta[j];
      }
    }
  }
  return SUCC;

}

//A: m by n matrix
//B: n by n matrix
//tC: l by n matrix
//D = A B^{-1} C^T
inline Int A_invB_transC(Mdoub A, pgpr_chol *chol_b, Mdoub C, Mdoub &D)
{
  Int m = A.nrows();
  Int n = A.ncols();
  Int l = C.nrows();
  Vdoub v(n);
  Vdoub beta(n);
  Vdoub Cn(n);
  D.assign(m, l, 0);
  for(Int i = 0; i < l; i++) {
    for(Int j = 0; j < n; j++) {
      Cn[j] = C[i][j];
    }
    chol_b->solve(Cn, beta);
    for(Int k = 0; k < m; k++) {
      for(Int j = 0; j < n; j++) {
        D[k][i] += A[k][j] * beta[j];
      }
    }
  }
  return SUCC;

}


//A: m by n matrix
//B: n by n matrix
//D = A B^{-1} A^T
inline Int A_invB_C(Mdoub A, pgpr_chol *chol_b, Mdoub &D)
{

  Int m = A.nrows();
  Int n = A.ncols();

  Vdoub v(n);
  Vdoub beta(n);
  Vdoub Cn(n);
  D.assign(m, m, 0);

  for(Int i = 0; i < m; i++) {

    for(Int j = 0; j < n; j++) {
      Cn[j] = A[i][j];
    }
    chol_b->elsolve(Cn, v);
    //D[i][i]
    for(Int j = 0; j < n; j++) {
      D[i][i] += v[j] * v[j];
    }
    //D[i][j] i!=j
    chol_b->solve(Cn, beta);
    for(Int k = i + 1; k < m; k++) {
      for(Int j = 0; j < n; j++) {
        D[k][i] += A[k][j] * beta[j];
      }
      D[i][k] = D[k][i];
    }
  }

  return SUCC;

}

//A: m by n matrix
//B: n by n matrix
//C: n by m matrix
//D = trace(A B^{-1} C)
inline Int trace_A_invB_C(Mdoub A, pgpr_chol *chol_b, Mdoub C, Vdoub &D)
{
  Int m = A.nrows();
  Int n = A.ncols();
  Int l = C.ncols();
  Vdoub v(n);
  Vdoub beta(n);
  Vdoub Cn(n);
  D.assign(m, 0);
  for(Int i = 0; i < l; i++) {
    //C(:,i)
    for(Int j = 0; j < n; j++) {
      Cn[j] = C[j][i];
    }
    //D[i]
    chol_b->solve(Cn, beta);
    for(Int j = 0; j < n; j++) {
      D[i] += A[i][j] * beta[j];
    }
  }
  return SUCC;

}

//A: m by n matrix
//B: n by n matrix
//D = trace(A B^{-1} A^T)
inline Int trace_A_invB_C(Mdoub A, pgpr_chol *chol_b, Vdoub &D)
{

  Int m = A.nrows();
  Int n = A.ncols();

  Vdoub v(n);
  Vdoub beta(n);
  Vdoub Cn(n);
  D.assign(m, 0);

  for(Int i = 0; i < m; i++) {

    for(Int j = 0; j < n; j++) {
      Cn[j] = A[i][j];
    }
    chol_b->elsolve(Cn, v);
    //D[i][i]
    for(Int j = 0; j < n; j++) {
      D[i] += v[j] * v[j];
    }
  }

  return SUCC;

}

//A: m by n matrix
//B: n by n matrix
//C: m by n matrix
//D = trace(A B^{-1} C^T)
inline Int trace_A_invB_transC(Mdoub A, pgpr_chol *chol_b, Mdoub C, Vdoub &D)
{
  Int m = A.nrows();
  Int n = A.ncols();
  Int l = C.nrows();
  if(m != l) {
    throw("something wrong with trace(AB^{-1}C^T)\n");
  }
  Vdoub v(n);
  Vdoub beta(n);
  Vdoub Cn(n);
  D.assign(m, 0);
  for(Int i = 0; i < m; i++) {
    //C(:,i)
    for(Int j = 0; j < n; j++) {
      Cn[j] = C[i][j];
    }
    //D[i]
    chol_b->solve(Cn, beta);
    for(Int j = 0; j < n; j++) {
      D[i] += A[i][j] * beta[j];
    }
  }
  return SUCC;

}

/**
 * @class pgpr_timer real-time timer
 * @brief This timer class can provide real-time measure (in seconds) incurred by a
 * block of running program 
 **/
class pgpr_timer
{
private:
  timespec st;
  timespec en;
  Doub elapsed;
public:
	///@brief constructor
  pgpr_timer() {
    elapsed = 0;
  }

	///@brief Start a timer
  inline void start() {
    // OS X does not have clock_gettime, use clock_get_time
#ifdef __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    st.tv_sec = mts.tv_sec;
    st.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, &st);
#endif
  }

	///@brief Stop a timer and compute the elapsed time in seconds 
	///@return elapsed time
  inline Doub end() {
#ifdef __MACH__
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    en.tv_sec = mts.tv_sec;
    en.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, &en);
#endif
    elapsed = (Doub (en.tv_sec - st.tv_sec))  + ((Doub(en.tv_nsec - st.tv_nsec)) / 1000000000);
    return elapsed;
  }
};

#endif
