/*
 * This example demonstrates the use of ISO_C_BINDING to call parpack (portability).
 * IMPORTANT: MPI communicators MUST be passed from C to Fortran using MPI_Comm_c2f.
 *
 * Just use arpack as you would have normally done, but, use *[ae]upd_c instead of *[ae]upd_.
 * The main advantage is that compiler checks (arguments) are performed at build time.
 * Note: to debug parpack, call debug_c.
 */

#include <iostream>
#include <cmath>
#include "mpi.h"
#include "parpack.hpp"
#include <complex.h> // creal, cimag.
#include <string>
#include "debug_c.hpp" // debug parpack.
#include "stat_c.hpp" // arpack statistics.
#include <stdio.h>

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * */
#ifndef BLASINT
#define BLASINT int
#endif

void sMatVec(float * x, float * y) {
  int i;
  for ( i = 0; i < 1000; ++i)
    y[i] = ((float) (i+1))*x[i];
};

int ss() {
  BLASINT ido = 0;
  std::string bmat("I");
  BLASINT N = 1000;
  std::string which("LM");
  BLASINT nev = 3;
  float tol = 0;
  float resid[N];
  BLASINT ncv = 2*nev+1;
  float V[ncv*N];
  BLASINT ldv = N;
  BLASINT iparam[11];
  BLASINT ipntr[14];
  float workd[3*N];
  bool rvec = true;
  std::string howmny("A");
  float* d = (float*) new float[(nev+1)];
  int select[ncv];
  float z[(N+1)*(nev+1)];
  BLASINT ldz = N+1;
  float sigma=0;
  int k;
  for (k=0; k < 3*N; ++k )
    workd[k] = 0;
  float workl[3*(ncv*ncv) + 6*ncv];
  for (k=0; k < 3*(ncv*ncv) + 6*ncv; ++k )
    workl[k] = 0;
  BLASINT lworkl = 3*(ncv*ncv) + 6*ncv;
  BLASINT info = 0;
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iparam[0] = 1;
  iparam[2] = 10*N;
  iparam[3] = 1;
  iparam[4] = 0; // number of ev found by arpack.
  iparam[6] = 1;

  MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);
  while(ido != 99) {
    /* call arpack like you would have, but, use ssaupd_c instead of ssaupd_ */
    pssaupd_c(MCW, ido, bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
              workd, workl, lworkl, info);

    sMatVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));
  }
  if (iparam[4] != nev) return 1; // check number of ev found by arpack.

  /* call arpack like you would have, but, use sseupd_c instead of sseupd_ */
  psseupd_c(MCW, rvec, howmny.c_str(), select, d, z, ldz, sigma,
            bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
            workd, workl, lworkl, info);
  int i;
  for (i = 0; i < nev; ++i) {
    std::cout << "rank " << rank << " - " << d[i] << std::endl;
    if(fabs(d[i] - (float)(1000-(nev-1)+i))>1e-1){
      delete [] d;
      return 1;
    }
  }
  delete [] d;
  return 0;
}

void cMatVec(float _Complex * x, float _Complex * y) {
  int i;
  for (i = 0; i < 1000; ++i)
    y[i] = x[i] * (i+1.0f + _Complex_I * (i+1.0f));
};

int cn() {
  BLASINT ido = 0;
  std::string bmat("I");
  BLASINT N = 1000;
  std::string which("LM");
  BLASINT nev = 1;
  float tol = 0;
  float _Complex resid[N];
  BLASINT ncv = 2*nev+1;
  float _Complex V[ncv*N];
  BLASINT ldv = N;
  BLASINT iparam[11];
  BLASINT ipntr[14];
  float _Complex workd[3*N];
  bool rvec = true;
  std::string howmny("A");
  float _Complex* d = (float _Complex*) new float _Complex[(nev+1)];
  int select[ncv];
  float _Complex z[(N+1)*(nev+1)];
  BLASINT ldz = N+1;
  float sigma=0;
  int k;
  for (k=0; k < 3*N; ++k )
    workd[k] = 0.f + I * 0.f;
  float _Complex workl[3*(ncv*ncv) + 6*ncv];
  for (k=0; k < 3*(ncv*ncv) + 6*ncv; ++k )
    workl[k] = 0.f + I * 0.f;
  BLASINT lworkl = 3*(ncv*ncv) + 6*ncv;
  float _Complex rwork[ncv];
  float _Complex workev[2*ncv];
  BLASINT info = 0;
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  iparam[0] = 1;
  iparam[2] = 10*N;
  iparam[3] = 1;
  iparam[4] = 0; // number of ev found by arpack.
  iparam[6] = 1;

  MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);
  while(ido != 99) {
    /* call arpack like you would have, but, use cnaupd_c instead of cnaupd_ */
    pcnaupd_c(MCW, ido, bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
              workd, workl, lworkl, rwork, info);

    cMatVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));
  }
  if (iparam[4] != nev) return 1; // check number of ev found by arpack.

  /* call arpack like you would have, but, use cneupd_c instead of cneupd_ */
  pcneupd_c(MCW, rvec, howmny.c_str(), select, d, z, ldz, sigma, workev,
            bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
            workd, workl, lworkl, rwork, info);
  int i;
  for (i = 0; i < nev; ++i) {
    std::cout << "rank " << rank << " - " << creal(d[i]) << " " << cimag(d[i]) << std::endl;
    if(fabs(creal(d[i]) - (float)(1000-i))>1e-1 || fabs(cimag(d[i]) - (float)(1000-i))>1e-1){
      delete [] d;
      return 1;
    }
  }
  delete [] d;
  return 0;
}

int main() {
  MPI_Init(NULL, NULL);

  sstats_c();
  int rc = ss(); // parpack without debug.
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rc != 0) return rc;
  int nopx_c, nbx_c, nrorth_c, nitref_c, nrstrt_c;
  float tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c;
  float tnaupd_c, tnaup2_c, tnaitr_c, tneigt_c, tngets_c, tnapps_c, tnconv_c;
  float tcaupd_c, tcaup2_c, tcaitr_c, tceigt_c, tcgets_c, tcapps_c, tcconv_c;
  float tmvopx_c, tmvbx_c, tgetv0_c, titref_c, trvec_c;
  stat_c(  nopx_c,    nbx_c, nrorth_c, nitref_c, nrstrt_c,
         tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c,
         tnaupd_c, tnaup2_c, tnaitr_c, tneigt_c, tngets_c, tnapps_c, tnconv_c,
         tcaupd_c, tcaup2_c, tcaitr_c, tceigt_c, tcgets_c, tcapps_c, tcconv_c,
         tmvopx_c,  tmvbx_c, tgetv0_c, titref_c,  trvec_c);
  std::cout << "Timers : nopx " << nopx_c << ", tmvopx " << tmvopx_c;
  std::cout <<        " - nbx " <<  nbx_c << ", tmvbx "  <<  tmvbx_c << std::endl;

  int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) std::cout << "------" << std::endl;

  debug_c(6, -6, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1); // set debug flags.
  rc = cn(); // parpack with debug.
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rc != 0) return rc;

  MPI_Finalize();
  return 0;
}
