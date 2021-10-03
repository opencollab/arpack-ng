#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "arpack.h"

// This test calls fortran from C the old-fashion cumbersome way.
// Note: icb_arpack_c tests the same kind of things using ICB.

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * We're using the non symmetric routines dnaupd and dneupd.
 * This is not efficient since the problem is
 * symmetric but is done to exhibit the bug.
 * */

void matVec(double *x, double *y) {
  int i;
  for (i = 0; i < 1000; ++i) y[i] = ((double)(i + 1)) * x[i];
};

int main() {
  a_int ido = 0;
  char bmat[] = "I";
  a_int N = 1000;
  char which[] = "LM";
  a_int nev = 9;
  double tol = 0;
  double resid[N];
  a_int ncv = 2 * nev + 1;
  double V[ncv * N];
  a_int ldv = N;
  a_int iparam[11];
  a_int ipntr[14];
  double workd[3 * N];
  a_int rvec = 1;
  char howmny[] = "A";
  double *dr = (double *)malloc((nev + 1) * sizeof(double));
  double *di = (double *)malloc((nev + 1) * sizeof(double));
  a_int select[3 * ncv];
  double z[(N + 1) * (nev + 1)];
  a_int ldz = N + 1;
  double sigmar = 0;
  double sigmai = 0;
  double workev[3 * ncv];
  int k;
  for (k = 0; k < 3 * N; ++k) workd[k] = 0;
  double workl[3 * (ncv * ncv) + 6 * ncv];
  for (k = 0; k < 3 * (ncv * ncv) + 6 * ncv; ++k) workl[k] = 0;
  a_int lworkl = 3 * (ncv * ncv) + 6 * ncv;
  a_int info = 0;

  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[6] = 1;

  dnaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
         workd, workl, lworkl, &info);

  while (ido == -1 || ido == 1) {
    matVec(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));

    dnaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam,
           ipntr, workd, workl, lworkl, &info);
  }

  dneupd_c(rvec, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, bmat,
         N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd,
         workl, lworkl, &info);
  int i;
  for (i = 0; i < nev; ++i) {
    printf("%f\n", dr[i]);
    if (fabs(dr[i] - (double)(1000 - i)) > 1e-6) {
      free(dr);
      free(di);
      exit(EXIT_FAILURE);
    }
  }
  free(dr);
  free(di);
  return 0;
}
