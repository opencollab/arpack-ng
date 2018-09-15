#ifndef __ARPACK_H__
#define __ARPACK_H__

#ifdef __cplusplus
extern "C" {
#endif

void cnaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, float _Complex* resid, int ncv, float _Complex* v, int ldv, int* iparam, int* ipntr, float _Complex* workd, float _Complex* workl, int lworkl, float* rwork, int* info);
void cneupd_c(bool rvec, char const* howmny, int const* select, float _Complex* d, float _Complex* z, int ldz, float _Complex sigma, float _Complex* workev, char const* bmat, int n, char const* which, int nev, float tol, float _Complex* resid, int ncv, float _Complex* v, int ldv, int* iparam, int* ipntr, float _Complex* workd, float _Complex* workl, int lworkl, float* rwork, int* info);
void dnaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void dneupd_c(bool rvec, char const* howmny, int const* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double * workev, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void dsaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void dseupd_c(bool rvec, char const* howmny, int const* select, double* d, double* z, int ldz, double sigma, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void snaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void sneupd_c(bool rvec, char const* howmny, int const* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float * workev, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void ssaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void sseupd_c(bool rvec, char const* howmny, int const* select, float* d, float* z, int ldz, float sigma, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void znaupd_c(int* ido, char const* bmat, int n, char const* which, int nev, double tol, double _Complex* resid, int ncv, double _Complex* v, int ldv, int* iparam, int* ipntr, double _Complex* workd, double _Complex* workl, int lworkl, double* rwork, int* info);
void zneupd_c(bool rvec, char const* howmny, int const* select, double _Complex* d, double _Complex* z, int ldz, double _Complex sigma, double _Complex* workev, char const* bmat, int n, char const* which, int nev, double tol, double _Complex* resid, int ncv, double _Complex* v, int ldv, int* iparam, int* ipntr, double _Complex* workd, double _Complex* workl, int lworkl, double* rwork, int* info);

#ifdef  __cplusplus
}
#endif

#endif
