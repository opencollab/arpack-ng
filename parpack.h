#ifndef __PARPACK_H__
#define __PARPACK_H__

/*
 * IMPORTANT: MPI communicators MUST be passed from C to Fortran using MPI_Comm_c2f.
 *            MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);
 */

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

void pcnaupd_c(MPI_Fint comm, int* ido, char const* bmat, int n, char const* which, int nev, float tol, float _Complex* resid, int ncv, float _Complex* v, int ldv, int* iparam, int* ipntr, float _Complex* workd, float _Complex* workl, int lworkl, float _Complex* rwork, int* info);
void pcneupd_c(MPI_Fint comm, bool rvec, char const* howmny, int const* select, float _Complex* d, float _Complex* z, int ldz, float _Complex sigma, float _Complex* workev, char const* bmat, int n, char const* which, int nev, float tol, float _Complex* resid, int ncv, float _Complex* v, int ldv, int* iparam, int* ipntr, float _Complex* workd, float _Complex* workl, int lworkl, float _Complex* rwork, int* info);
void pdnaupd_c(MPI_Fint comm, int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void pdneupd_c(MPI_Fint comm, bool rvec, char const* howmny, int const* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void pdsaupd_c(MPI_Fint comm, int* ido, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void pdseupd_c(MPI_Fint comm, bool rvec, char const* howmny, int const* select, double* d, double* z, int ldz, double sigma, char const* bmat, int n, char const* which, int nev, double tol, double* resid, int ncv, double* v, int ldv, int* iparam, int* ipntr, double* workd, double* workl, int lworkl, int* info);
void psnaupd_c(MPI_Fint comm, int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void psneupd_c(MPI_Fint comm, bool rvec, char const* howmny, int const* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void pssaupd_c(MPI_Fint comm, int* ido, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void psseupd_c(MPI_Fint comm, bool rvec, char const* howmny, int const* select, float* d, float* z, int ldz, float sigma, char const* bmat, int n, char const* which, int nev, float tol, float* resid, int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd, float* workl, int lworkl, int* info);
void pznaupd_c(MPI_Fint comm, int* ido, char const* bmat, int n, char const* which, int nev, double tol, double _Complex* resid, int ncv, double _Complex* v, int ldv, int* iparam, int* ipntr, double _Complex* workd, double _Complex* workl, int lworkl, double _Complex* rwork, int* info);
void pzneupd_c(MPI_Fint comm, bool rvec, char const* howmny, int const* select, double _Complex* d, double _Complex* z, int ldz, double _Complex sigma, double _Complex* workev, char const* bmat, int n, char const* which, int nev, double tol, double _Complex* resid, int ncv, double _Complex* v, int ldv, int* iparam, int* ipntr, double _Complex* workd, double _Complex* workl, int lworkl, double _Complex* rwork, int* info);

#ifdef  __cplusplus
}
#endif

#endif
