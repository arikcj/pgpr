/**
 * @file pgpr_type.h
 * @brief This file provides important macros, templates, basic data
 * types (e.g., vector, matrix). 
 * 
 * @author CHEN jie, arik.cj@gmail.com
 * 
 * @version 1.0
 **/

#ifndef _PGPR_TYPE_H_
#define _PGPR_TYPE_H_
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <cfloat>
#include <stdarg.h>

using namespace std;
//macro
#define SUCC 0
#define FAIL -1
#define FALSE 0
#define TRUE 1


///@class pgpr_vector Vector class
///@brief Vector class 
template <class T>
class pgpr_vector
{
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	pgpr_vector();
	explicit pgpr_vector(int n);		// Zero-based array
	pgpr_vector(int n, const T &a);	//initialize to constant value
	pgpr_vector(int n, const T *a);	// Initialize to array
	pgpr_vector(const pgpr_vector &rhs);	// Copy constructor
	pgpr_vector & operator=(const pgpr_vector &rhs);	//assignment
	pgpr_vector & operator-(const pgpr_vector &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	~pgpr_vector();
};

// pgpr_vector definitions

template <class T>
pgpr_vector<T>::pgpr_vector() : nn(0), v(NULL) {}

template <class T>
pgpr_vector<T>::pgpr_vector(int n) : nn(n), v(n > 0 ? new T[n] : NULL) {}

template <class T>
pgpr_vector<T>::pgpr_vector(int n, const T& a) : nn(n), v(n > 0 ? new T[n] : NULL)
{
	for(int i = 0; i < n; i++) v[i] = a;
}

template <class T>
pgpr_vector<T>::pgpr_vector(int n, const T *a) : nn(n), v(n > 0 ? new T[n] : NULL)
{
	for(int i = 0; i < n; i++) v[i] = *a++;
}

template <class T>
pgpr_vector<T>::pgpr_vector(const pgpr_vector<T> &rhs) : nn(rhs.nn), v(nn > 0 ? new T[nn] : NULL)
{
	for(int i = 0; i < nn; i++) v[i] = rhs[i];
}

template <class T>
pgpr_vector<T> & pgpr_vector<T>::operator=(const pgpr_vector<T> &rhs){
/* postcondition: normal assignment via copying has been performed;
**if vector and rhs were different sizes, vector
**		has been resized to match the size of rhs
*/
	if (this != &rhs) {
		if (nn != rhs.nn) {
			if (v != NULL) delete [] (v);
			nn = rhs.nn;
			v = nn > 0 ? new T[nn] : NULL;
		}
		for (int i = 0; i < nn; i++)
			v[i] = rhs[i];
	}
	return *this;
}
template <class T>
inline T & pgpr_vector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw("pgpr_vector subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline const T & pgpr_vector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw("pgpr_vector subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline int pgpr_vector<T>::size() const
{
	return nn;
}

template <class T>
void pgpr_vector<T>::resize(int newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

template <class T>
void pgpr_vector<T>::assign(int newn, const T& a)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i = 0; i < nn; i++) v[i] = a;
}

template <class T>
pgpr_vector<T>::~pgpr_vector()
{
	if (v != NULL) delete[] (v);
}

// end of pgpr_vector definitions

///@class pgpr_matrix Matrix class
///@brief Matrix class 
template <class T>
class pgpr_matrix
{
private:
	int nn;
	int mm;
	T **v;
public:
	pgpr_matrix();
	pgpr_matrix(int n, int m);			// Zero-based array
	pgpr_matrix(int n, int m, const T &a);	//Initialize to constant
	pgpr_matrix(int n, int m, const T *a);	// Initialize to array
	pgpr_matrix(const pgpr_matrix &rhs);		// Copy constructor
	pgpr_matrix & operator=(const pgpr_matrix &rhs);	//assignment
	pgpr_matrix & operator+=(const pgpr_matrix &rhs);	//self addition
	pgpr_matrix & operator-=(const pgpr_matrix &rhs);	//self subtraction
	const pgpr_matrix operator+(const pgpr_matrix &b); //reload matrix addition
	const pgpr_matrix operator-(const pgpr_matrix &b); //reload matrix addition
	const pgpr_matrix operator*(const pgpr_matrix &b); //reload matrix subtraction
	const pgpr_matrix operator~(); //define: A~ = A^-1
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	~pgpr_matrix();
};

template <class T>
pgpr_matrix<T>::pgpr_matrix() : nn(0), mm(0), v(NULL) {}

template <class T>
pgpr_matrix<T>::pgpr_matrix(int n, int m) : nn(n), mm(m), v(n > 0 ? new T*[n] : NULL)
{
	int i, nel = m * n;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
}

template <class T>
pgpr_matrix<T>::pgpr_matrix(int n, int m, const T &a) : nn(n), mm(m), v(n > 0 ? new T*[n] : NULL)
{
	int i, j, nel = m * n;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
	for (i = 0; i < n; i++) for (j = 0; j < m; j++) v[i][j] = a;
}

template <class T>
pgpr_matrix<T>::pgpr_matrix(int n, int m, const T *a) : nn(n), mm(m), v(n > 0 ? new T*[n] : NULL)
{
	int i, j, nel = m * n;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < n; i++) v[i] = v[i - 1] + m;
	for (i = 0; i < n; i++) for (j = 0; j < m; j++) v[i][j] = *a++;
}

template <class T>
pgpr_matrix<T>::pgpr_matrix(const pgpr_matrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn > 0 ? new T*[nn] : NULL)
{
	int i, j, nel = mm * nn;
	if (v) v[0] = nel > 0 ? new T[nel] : NULL;
	for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
pgpr_matrix<T> & pgpr_matrix<T>::operator=(const pgpr_matrix<T> &rhs)
/* postcondition: normal assignment via copying has been performed;
**		if matrix and rhs were different sizes, matrix
**		has been resized to match the size of rhs
*/
{
	if (this != &rhs) {
		int i, j, nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn = rhs.nn;
			mm = rhs.mm;
			v = nn > 0 ? new T*[nn] : NULL;
			nel = mm * nn;
			if (v) v[0] = nel > 0 ? new T[nel] : NULL;
			for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
		}
		for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}
template <class T>
pgpr_matrix<T> & pgpr_matrix<T>::operator+=(const pgpr_matrix<T> &rhs){
	int i, j, nel;
	if (nn != rhs.nn){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}
	if (mm != rhs.mm){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] += rhs[i][j];
	return *this;
}

template <class T>
pgpr_matrix<T> & pgpr_matrix<T>::operator-=(const pgpr_matrix<T> &rhs)
{
	int i, j, nel;
	if (nn != rhs.nn){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}
	if (mm != rhs.mm){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] -= rhs[i][j];
	return *this;
}

template <class T>
inline T* pgpr_matrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw("pgpr_matrix subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline const T* pgpr_matrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= nn) {
		throw("pgpr_matrix subscript out of bounds");
	}
#endif
	return v[i];
}

template <class T>
inline const pgpr_matrix<T>  pgpr_matrix<T>::operator+(const pgpr_matrix<T> &b){
	if (nn != b.nrows()){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}

	if (mm != b.ncols()){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}
	pgpr_matrix<T> sum(nn, mm);
	int i, j;
	for(i = 0; i < nn; i ++)
		for(j = 0; j < mm; j ++)
			sum[i][j] = v[i][j] + b [i][j];
	return sum;		
}


template <class T>
inline const pgpr_matrix<T>  pgpr_matrix<T>::operator*(const pgpr_matrix<T> &b){
	if (mm != b.nn){
		throw("pgpr_matrix multiplication,dimension not compatible.");
	}
	pgpr_matrix<T> sum(nn, mm);
	int i, j;
	for(i = 0; i < nn; i ++)
		for(j = 0; j < mm; j ++)
			sum[i][j] = v[i][j] - b [i][j];
	return sum;		
}

template <class T>
inline const pgpr_matrix<T>  pgpr_matrix<T>::operator-(const pgpr_matrix<T> &b){
	if (nn != b.nrows()){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}

	if (mm != b.ncols()){
		throw("pgpr_matrix addition, not consistant in row dimension");
	}
	pgpr_matrix<T> sum(nn, mm);
	int i, j;
	for(i = 0; i < nn; i ++)
		for(j = 0; j < mm; j ++)
			sum[i][j] = v[i][j] - b [i][j];
	return sum;		
}

template <class T>
inline int pgpr_matrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline int pgpr_matrix<T>::ncols() const
{
	return mm;
}

template <class T>
void pgpr_matrix<T>::resize(int newn, int newm)
{
	int i, nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn > 0 ? new T*[nn] : NULL;
		nel = mm * nn;
		if (v) v[0] = nel > 0 ? new T[nel] : NULL;
		for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
	}
}

template <class T>
void pgpr_matrix<T>::assign(int newn, int newm, const T& a)
{
	int i, j, nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn > 0 ? new T*[nn] : NULL;
		nel = mm * nn;
		if (v) v[0] = nel > 0 ? new T[nel] : NULL;
		for (i = 1; i < nn; i++) v[i] = v[i - 1] + mm;
	}
	for (i = 0; i < nn; i++) for (j = 0; j < mm; j++) v[i][j] = a;
}

template <class T>
pgpr_matrix<T>::~pgpr_matrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}


// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const pgpr_vector<Bool> Vbool_I;
typedef pgpr_vector<Bool> Vbool, Vbool_O, Vbool_IO;

typedef const pgpr_vector<Int> Vint_I;
typedef pgpr_vector<Int> Vint, Vint_O, Vint_IO;

typedef const pgpr_vector<Uint> Vuint_I;
typedef pgpr_vector<Uint> Vuint, Vuint_O, Vuint_IO;

typedef const pgpr_vector<Doub> Vdoub_I;
typedef pgpr_vector<Doub> Vdoub, Vdoub_O, Vdoub_IO;

// matrix types

typedef const pgpr_matrix<Int> Mint_I;
typedef pgpr_matrix<Int> Mint, Mint_O, Mint_IO;

typedef const pgpr_matrix<Uint> Muint_I;
typedef pgpr_matrix<Uint> Muint, Muint_O, Muint_IO;

typedef const pgpr_matrix<Doub> Mdoub_I;
typedef pgpr_matrix<Doub> Mdoub, Mdoub_O, Mdoub_IO;

#endif /* _PGPR_H_ */

