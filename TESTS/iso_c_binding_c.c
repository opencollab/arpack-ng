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

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * */
#ifndef BLASINT
#define BLASINT int
#endif

extern void dsaupd_c(BLASINT * ido, char * bmat, BLASINT n, char * which, BLASINT nev,
                     double tol, double * resid, BLASINT ncv, double * v,
                     BLASINT ldv, BLASINT * iparam, BLASINT * ipntr, double * workd,
                     double * workl, BLASINT lworkl, BLASINT * info);

extern void dseupd_c(bool rvec, char * howmny, bool * select, double * d, double * z, int ldz, double sigma,
                     char * bmat, BLASINT n, char * which, BLASINT nev,
                     double tol, double * resid, BLASINT ncv, double * v,
                     BLASINT ldv, BLASINT * iparam, BLASINT * ipntr, double * workd,
                     double * workl, BLASINT lworkl, BLASINT * info);

void matVec(double * x, double * y) {
  int i;
  for ( i = 0; i < 1000; ++i)
    y[i] = ((double) (i+1))*x[i];
};

int main() {
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
  iparam[6] = 1;

  /* call arpack like you would have, but, use dsaupd_c instead of dsaupd_ */
  dsaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);

  while(ido == 1) {

    matVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));

    /* call arpack like you would have, but, use dsaupd_c instead of dsaupd_ */
    dsaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, &info);
  }

  /* call arpack like you would have, but, use dseupd_c instead of dseupd_ */
  dseupd_c(rvec, howmny, select, d, z, ldz, sigma,
           bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);
  int i;
  for (i = 0; i < nev; ++i) {
    printf("%f\n", d[i]);
    if(fabs(d[i] - (double)(1000-(nev-1)+i))>1e-6){
      free(d);
      exit(EXIT_FAILURE);
    }
  }
  free(d);
  return 0;
}
