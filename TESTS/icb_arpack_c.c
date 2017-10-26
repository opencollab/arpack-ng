/*
 * This example demonstrates the use of ISO_C_BINDING to call arpack (portability).
 *
 * Just use arpack as you would have normally done, but, use *[ae]upd_c instead of *[ae]upd_.
 * The main advantage is that compiler checks (arguments) are performed at build time.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h> // bool.
#include "arpack.h"
#include <complex.h> // creal, cimag.

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * */
#ifndef BLASINT
#define BLASINT int
#endif

void dMatVec(double * x, double * y) {
  int i;
  for ( i = 0; i < 1000; ++i)
    y[i] = ((double) (i+1))*x[i];
};

int ds() {
  BLASINT ido = 0;
  char bmat[] = "I";
  BLASINT N = 1000;
  char which[] = "LM";
  BLASINT nev = 9;
  double tol = 0;
  double resid[N];
  BLASINT ncv = 2*nev+1;
  double V[ncv*N];
  BLASINT ldv = N;
  BLASINT iparam[11];
  BLASINT ipntr[14];
  double workd[3*N];
  bool rvec = true;
  char howmny[] = "A";
  double* d = (double*) malloc((nev+1)*sizeof(double));
  bool select[3*ncv];
  double z[(N+1)*(nev+1)];
  BLASINT ldz = N+1;
  double sigma=0;
  int k;
  for (k=0; k < 3*N; ++k )
    workd[k] = 0;
  double workl[3*(ncv*ncv) + 6*ncv];
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
    /* call arpack like you would have, but, use dsaupd_c instead of dsaupd_ */
    dsaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, &info);

    dMatVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));
  }
  if (iparam[4] != nev) return 1; // check number of ev found by arpack.

  /* call arpack like you would have, but, use dseupd_c instead of dseupd_ */
  dseupd_c(rvec, howmny, select, d, z, ldz, sigma,
           bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);
  int i;
  for (i = 0; i < nev; ++i) {
    printf("%f\n", d[i]);
    if(fabs(d[i] - (double)(1000-(nev-1)+i))>1e-6){
      free(d);
      return 1;
    }
  }
  free(d);
  return 0;
}

void zMatVec(double _Complex * x, double _Complex * y) {
  int i;
  for (i = 0; i < 1000; ++i)
    y[i] = (i+1.)*creal(x[i]) + I * (i+1.)*cimag(x[i]);
};

int zn() {
  BLASINT ido = 0;
  char bmat[] = "I";
  BLASINT N = 1000;
  char which[] = "LM";
  BLASINT nev = 9;
  double tol = 0;
  double _Complex resid[N];
  BLASINT ncv = 2*nev+1;
  double _Complex V[ncv*N];
  BLASINT ldv = N;
  BLASINT iparam[11];
  BLASINT ipntr[14];
  double _Complex workd[3*N];
  bool rvec = true;
  char howmny[] = "A";
  double _Complex* d = (double _Complex*) malloc((nev+1)*sizeof(double _Complex));
  bool select[3*ncv];
  double _Complex z[(N+1)*(nev+1)];
  BLASINT ldz = N+1;
  double sigma=0;
  int k;
  for (k=0; k < 3*N; ++k )
    workd[k] = 0;
  double _Complex workl[3*(ncv*ncv) + 6*ncv];
  for (k=0; k < 3*(ncv*ncv) + 6*ncv; ++k )
    workl[k] = 0;
  BLASINT lworkl = 3*(ncv*ncv) + 6*ncv;
  double _Complex rwork[ncv];
  double _Complex workev[2*ncv];
  BLASINT info = 0;

  iparam[0] = 1;
  iparam[2] = 10*N;
  iparam[3] = 1;
  iparam[4] = 0; // number of ev found by arpack.
  iparam[6] = 1;

  while(ido != 99) {
    /* call arpack like you would have, but, use znaupd_c instead of znaupd_ */
    znaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, rwork, &info);

    zMatVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));
  }
  if (iparam[4] != nev) return 1; // check number of ev found by arpack.

  /* call arpack like you would have, but, use zneupd_c instead of zneupd_ */
  zneupd_c(rvec, howmny, select, d, z, ldz, sigma, workev,
           bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, rwork, &info);
  int i;
  for (i = 0; i < nev; ++i) {
    printf("%f %f\n", creal(d[i]), cimag(d[i]));
    if(fabs(creal(d[i]) - (double)(1000-i))>1e-6){
      free(d);
      return 1;
    }
  }
  free(d);
  return 0;
}

int main() {
  if (ds() != 0) return 1;
  printf("------\n");
  return zn();
}
