/*
 * This example demonstrates the use of ISO_C_BINDING to call arpack (portability).
 *
 * Just use arpack as you would have normally done, but, use *[ae]upd_c instead of *[ae]upd_.
 * The main advantage is that compiler checks (arguments) are performed at build time.
 */

#include <iostream>
#include <cmath>
#include "arpack.hpp"

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * */
#ifndef BLASINT
#define BLASINT int
#endif

void matVec(float * x, float * y) {
  int i;
  for ( i = 0; i < 1000; ++i)
    y[i] = ((float) (i+1))*x[i];
};

int main() {
  BLASINT ido = 0;
  char bmat[] = "I";
  BLASINT N = 1000;
  char which[] = "LM";
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
  char howmny[] = "A";
  float* d = (float*) new float[(nev+1)];
  bool select[3*ncv];
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
  iparam[6] = 1;

  /* call arpack like you would have, but, use ssaupd_c instead of ssaupd_ */
  ssaupd_c(ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, info);

  while(ido == 1) {

    matVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));

    /* call arpack like you would have, but, use ssaupd_c instead of ssaupd_ */
    ssaupd_c(ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, info);
  }

  /* call arpack like you would have, but, use sseupd_c instead of sseupd_ */
  sseupd_c(rvec, howmny, select, d, z, ldz, sigma,
           bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
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
