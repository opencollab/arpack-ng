#ifndef __ARPACK_HPP__
#define __ARPACK_HPP__

/*
 * From C++, arpack does not exist.
 * Arpack is Fortran. ISO_C_BINDING is a gateway from Fortran to C, not C++. But, C++ can "get back" to C.
 * This is why you find C types in the arpack.hpp (C++ header) to be as safe as possible.
 */

extern "C" void ssaupd_c(int & ido, const char * bmat, int n, const char * which, int nev,
                         float tol, float * resid, int ncv, float * v,
                         int ldv, int * iparam, int * ipntr, float * workd,
                         float * workl, int lworkl, int & info);

extern "C" void sseupd_c(bool rvec, const char * howmny, int * select, float * d, float * z, int ldz, float sigma,
                         const char * bmat, int n, const char * which, int nev,
                         float tol, float * resid, int ncv, float * v,
                         int ldv, int * iparam, int * ipntr, float * workd,
                         float * workl, int lworkl, int & info);

extern "C" void dsaupd_c(int & ido, const char * bmat, int n, const char * which, int nev,
                         double tol, double * resid, int ncv, double * v,
                         int ldv, int * iparam, int * ipntr, double * workd,
                         double * workl, int lworkl, int & info);

extern "C" void dseupd_c(bool rvec, const char * howmny, int * select, double * d, double * z, int ldz, double sigma,
                         const char * bmat, int n, const char * which, int nev,
                         double tol, double * resid, int ncv, double * v,
                         int ldv, int * iparam, int * ipntr, double * workd,
                         double * workl, int lworkl, int & info);

extern "C" void snaupd_c(int & ido, const char * bmat, int n, const char * which, int nev,
                         float tol, float * resid, int ncv, float * v,
                         int ldv, int * iparam, int * ipntr, float * workd,
                         float * workl, int lworkl, int & info);

extern "C" void sneupd_c(bool rvec, const char * howmny, int * select, float * dr, float * di, float * z, int ldz, float sigmar, float sigmai,
                         const char * bmat, int n, const char * which, int nev,
                         float tol, float * resid, int ncv, float * v,
                         int ldv, int * iparam, int * ipntr, float * workd,
                         float * workl, int lworkl, int & info);

extern "C" void dnaupd_c(int & ido, const char * bmat, int n, const char * which, int nev,
                         double tol, double * resid, int ncv, double * v,
                         int ldv, int * iparam, int * ipntr, double * workd,
                         double * workl, int lworkl, int & info);

extern "C" void dneupd_c(bool rvec, const char * howmny, int * select, double * dr, double * di, double * z, int ldz, double sigmar, double sigmai,
                         const char * bmat, int n, const char * which, int nev,
                         double tol, double * resid, int ncv, double * v,
                         int ldv, int * iparam, int * ipntr, double * workd,
                         double * workl, int lworkl, int & info);

extern "C" void cnaupd_c(int & ido, const char * bmat, int n, const char * which, int nev,
                         float tol, float _Complex * resid, int ncv, float _Complex * v,
                         int ldv, int * iparam, int * ipntr, float _Complex * workd,
                         float _Complex * workl, int lworkl, float _Complex * rwork, int & info);

extern "C" void cneupd_c(bool rvec, const char * howmny, int * select,
                         float _Complex * d, float _Complex * z, int ldz, float _Complex sigma, float _Complex * workev,
                         const char * bmat, int n, const char * which, int nev,
                         float tol, float _Complex * resid, int ncv, float _Complex * v,
                         int ldv, int * iparam, int * ipntr, float _Complex * workd,
                         float _Complex * workl, int lworkl, float _Complex * rwork, int & info);

extern "C" void znaupd_c(int & ido, const char * bmat, int n, const char * which, int nev,
                         double tol, double _Complex * resid, int ncv, double _Complex * v,
                         int ldv, int * iparam, int * ipntr, double _Complex * workd,
                         double _Complex * workl, int lworkl, double _Complex * rwork, int & info);

extern "C" void zneupd_c(bool rvec, const char * howmny, int * select,
                         double _Complex * d, double _Complex * z, int ldz, double _Complex sigma, double _Complex * workev,
                         const char * bmat, int n, const char * which, int nev,
                         double tol, double _Complex * resid, int ncv, double _Complex * v,
                         int ldv, int * iparam, int * ipntr, double _Complex * workd,
                         double _Complex * workl, int lworkl, double _Complex * rwork, int & info);

#endif
