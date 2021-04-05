/*
 * This example demonstrates the use of ISO_C_BINDING to call arpack
 * (portability).
 *
 * Just use arpack as you would have normally done, but, use *[ae]upd_c instead
 * of *[ae]upd_. The main advantage is that compiler checks (arguments) are
 * performed at build time. Note: to debug arpack, call debug_c.
 */

#include <complex.h>  // creal, cimag.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "arpack.h"
#include "debug_c.h"  // debug arpack.
#include "stat_c.h"   // arpack statistics.

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * */

void dMatVec(double* x, double* y) {
  int i;
  for (i = 0; i < 1000; ++i) y[i] = ((double)(i + 1)) * x[i];
};

int ds() {
  a_int ido = 0;
  char bmat[] = "I";
  a_int N = 1000;
  char which[] = "LM";
  a_int nev = 9;
  double tol = 0.000001; // small tol => more stable checks after EV computation.
  double resid[N];
  a_int ncv = 2 * nev + 1;
  double V[ncv * N];
  a_int ldv = N;
  a_int iparam[11];
  a_int ipntr[14];
  double workd[3 * N];
  a_int rvec = 1;
  char howmny[] = "A";
  double* d = (double*)malloc((nev + 1) * sizeof(double));
  a_int select[ncv];
  int i; // C99 compliant.
  for (i = 0; i < ncv; i++) select[i] = 1;
  double z[(N + 1) * (nev + 1)];
  a_int ldz = N + 1;
  double sigma = 0;
  int k;
  for (k = 0; k < 3 * N; ++k) workd[k] = 0;
  double workl[3 * (ncv * ncv) + 6 * ncv];
  for (k = 0; k < 3 * (ncv * ncv) + 6 * ncv; ++k) workl[k] = 0;
  a_int lworkl = 3 * (ncv * ncv) + 6 * ncv;
  a_int info = 0;

  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  while (ido != 99) {
    /* call arpack like you would have, but, use dsaupd_c instead of dsaupd_ */
    dsaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, &info);

    dMatVec(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  }
  if (iparam[4] != nev) {printf("Error: iparam[4] %d, nev %d\n", iparam[4], nev); return 1;} // check number of ev found by arpack.

  /* call arpack like you would have, but, use dseupd_c instead of dseupd_ */
  dseupd_c(rvec, howmny, select, d, z, ldz, sigma, bmat, N, which, nev, tol,
           resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, &info);
  for (i = 0; i < nev; ++i) {
    double val = d[i];
    double ref = (N-(nev-1)+i);
    double eps = fabs(val - ref);
    printf("%f - %f - %f\n", val, ref, eps);

    /*eigen value order: smallest -> biggest*/
    if(eps>1.e-05){
      free(d);
      return 1;
    }
  }
  free(d);
  return 0;
}

void zMatVec(double _Complex* x, double _Complex* y) {
  int i;
  for (i = 0; i < 1000; ++i) y[i] = x[i] * (i + 1.0 + _Complex_I * (i + 1.0));
};

int zn() {
  a_int ido = 0;
  char bmat[] = "I";
  a_int N = 1000;
  char which[] = "LM";
  a_int nev = 9;
  double tol = 0.000001; // small tol => more stable checks after EV computation.
  double _Complex resid[N];
  a_int ncv = 2 * nev + 1;
  double _Complex V[ncv * N];
  a_int ldv = N;
  a_int iparam[11];
  a_int ipntr[14];
  double _Complex workd[3 * N];
  a_int rvec = 0;
  char howmny[] = "A";
  double _Complex* d =
      (double _Complex*)malloc((nev + 1) * sizeof(double _Complex));
  a_int select[ncv];
  int i; // C99 compliant.
  for (i = 0; i < ncv; i++) select[i] = 1;
  double _Complex z[(N + 1) * (nev + 1)];
  a_int ldz = N + 1;
  double _Complex sigma = 0. + I * 0.;
  int k;
  for (k = 0; k < 3 * N; ++k) workd[k] = 0;
  double _Complex workl[3 * (ncv * ncv) + 6 * ncv];
  for (k = 0; k < 3 * (ncv * ncv) + 6 * ncv; ++k) workl[k] = 0;
  a_int lworkl = 3 * (ncv * ncv) + 6 * ncv;
  double rwork[ncv];
  double _Complex workev[2 * ncv];
  a_int info = 0;

  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  while (ido != 99) {
    /* call arpack like you would have, but, use znaupd_c instead of znaupd_ */
    znaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, rwork, &info);

    zMatVec(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  }
  if (iparam[4] != nev) {printf("Error: iparam[4] %d, nev %d\n", iparam[4], nev); return 1;} // check number of ev found by arpack.

  /* call arpack like you would have, but, use zneupd_c instead of zneupd_ */
  zneupd_c(rvec, howmny, select, d, z, ldz, sigma, workev, bmat, N, which, nev,
           tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, rwork,
           &info);
  for (i = 0; i < nev; ++i) {
    double rval = creal(d[i]);
    double rref = (N-(nev-1)+i);
    double reps = fabs(rval - rref);
    double ival = cimag(d[i]);
    double iref = (N-(nev-1)+i);
    double ieps = fabs(ival - iref);
    printf("%f %f - %f %f - %f %f\n", rval, ival, rref, iref, reps, ieps);

    /*eigen value order: smallest -> biggest*/
    if(reps>1.e-05 || ieps>1.e-05){
      free(d);
      return 1;
    }
  }
  free(d);
  return 0;
}

int main() {
  sstats_c();
  int rc = ds();  // arpack without debug.
  if (rc != 0) return rc;
  a_int nopx_c, nbx_c, nrorth_c, nitref_c, nrstrt_c;
  float tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c;
  float tnaupd_c, tnaup2_c, tnaitr_c, tneigt_c, tngets_c, tnapps_c, tnconv_c;
  float tcaupd_c, tcaup2_c, tcaitr_c, tceigt_c, tcgets_c, tcapps_c, tcconv_c;
  float tmvopx_c, tmvbx_c, tgetv0_c, titref_c, trvec_c;
  stat_c(&nopx_c, &nbx_c, &nrorth_c, &nitref_c, &nrstrt_c, &tsaupd_c, &tsaup2_c,
         &tsaitr_c, &tseigt_c, &tsgets_c, &tsapps_c, &tsconv_c, &tnaupd_c,
         &tnaup2_c, &tnaitr_c, &tneigt_c, &tngets_c, &tnapps_c, &tnconv_c,
         &tcaupd_c, &tcaup2_c, &tcaitr_c, &tceigt_c, &tcgets_c, &tcapps_c,
         &tcconv_c, &tmvopx_c, &tmvbx_c, &tgetv0_c, &titref_c, &trvec_c);
  printf("Timers : nopx %d, tmvopx %f - nbx %d, tmvbx %f\n", nopx_c, tmvopx_c,
         nbx_c, tmvbx_c);

  printf("------\n");

  // clang-format off
  debug_c(6, -6, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1); // set debug flags.
  // clang-format on
  rc = zn();  // arpack with debug.

  return rc;
}
