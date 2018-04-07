/*
 * This example demonstrates the use of ISO_C_BINDING to call arpack (portability).
 *
 * Just use arpack as you would have normally done, but, use *[ae]upd_c instead of *[ae]upd_.
 * The main advantage is that compiler checks (arguments) are performed at build time.
 * Note: to debug arpack, call debug_c.
 */

#include <iostream>
#include <cmath>
#include "arpack.hpp"
#include <complex.h> // creal, cimag.
#include <string>
#include "debug_c.hpp" // debug arpack.

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
  BLASINT nev = 9;
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

  iparam[0] = 1;
  iparam[2] = 10*N;
  iparam[3] = 1;
  iparam[4] = 0; // number of ev found by arpack.
  iparam[6] = 1;

  while(ido != 99) {
    /* call arpack like you would have, but, use ssaupd_c instead of ssaupd_ */
    ssaupd_c(ido, bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, info);

    sMatVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));
  }
  if (iparam[4] != nev) return 1; // check number of ev found by arpack.

  /* call arpack like you would have, but, use sseupd_c instead of sseupd_ */
  sseupd_c(rvec, howmny.c_str(), select, d, z, ldz, sigma,
           bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, info);
  int i;
  for (i = 0; i < nev; ++i) {
    std::cout << d[i] << std::endl;
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
  BLASINT nev = 9;
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
    workd[k] = 0;
  float _Complex workl[3*(ncv*ncv) + 6*ncv];
  for (k=0; k < 3*(ncv*ncv) + 6*ncv; ++k )
    workl[k] = 0;
  BLASINT lworkl = 3*(ncv*ncv) + 6*ncv;
  float _Complex rwork[ncv];
  float _Complex workev[2*ncv];
  BLASINT info = 0;

  iparam[0] = 1;
  iparam[2] = 10*N;
  iparam[3] = 1;
  iparam[4] = 0; // number of ev found by arpack.
  iparam[6] = 1;

  while(ido != 99) {
    /* call arpack like you would have, but, use cnaupd_c instead of cnaupd_ */
    cnaupd_c(ido, bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, rwork, info);

    cMatVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));
  }
  if (iparam[4] != nev) return 1; // check number of ev found by arpack.

  /* call arpack like you would have, but, use cneupd_c instead of cneupd_ */
  cneupd_c(rvec, howmny.c_str(), select, d, z, ldz, sigma, workev,
           bmat.c_str(), N, which.c_str(), nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, rwork, info);
  int i;
  for (i = 0; i < nev; ++i) {
    std::cout << creal(d[i]) << " " << cimag(d[i]) << std::endl;
    if(fabs(creal(d[i]) - (float)(1000-i))>1e-1 || fabs(cimag(d[i]) - (float)(1000-i))>1e-1){
      delete [] d;
      return 1;
    }
  }
  delete [] d;
  return 0;
}

int main() {
  int rc = ss(); // arpack without debug.
  if (rc != 0) return rc;

  std::cout << "------" << std::endl;

  debug_c(6, -6, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1); // set debug flags.
  rc = cn(); // arpack with debug.

  return rc;
}
