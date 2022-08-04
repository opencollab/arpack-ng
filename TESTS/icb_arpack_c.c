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

void dMatVec(const double* x, double* y) {
  int i;
  for (i = 0; i < 1000; ++i) y[i] = ((double)(i + 1)) * x[i];
};

int ds() {
  const a_int N      = 1000;
  const a_int nev    = 9;
  const a_int ncv    = 2 * nev + 1;
  const a_int ldv    = N;
  const a_int ldz    = N;
  const a_int lworkl = ncv * (ncv + 8);
  const a_int rvec   = 1;      // need eigenvectors

  const double tol = 0.000001; // small tol => more stable checks after EV computation.
  const double sigma = 0;      // not referenced in this mode
  
  double resid[N];
  double V[ldv * ncv];
  double z[ldz * nev];
  double d[nev];
  double workd[3 * N];
  double workl[lworkl];
  a_int select[ncv]; // since HOWMNY = 'A', only used as workspace here

  a_int iparam[11], ipntr[11];
  iparam[0] = 1;       // ishift
  iparam[2] = 10 * N;  // on input: maxit; on output: actual iteration
  iparam[3] = 1;       // NB, only 1 allowed
  iparam[6] = 1;       // mode

  char bmat[]   = "I";
  char which[]  = "LM";
  char howmny[] = "A";

  a_int info = 0, ido = 0;
  do {
    dsaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, &info);

    dMatVec(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  } while (ido == 1 || ido == -1);

  // check info and number of ev found by arpack.
  if (info < 0 || iparam[4] < nev) {
    printf("Error in saupd: iparam[4] %d, nev %d, info %d\n", iparam[4], nev, info);
    return 1;
  }

  dseupd_c(rvec, howmny, select, d, z, ldz, sigma, bmat, N, which, nev, tol,
           resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, &info);
  if (info < 0) {
    printf("Error in seupd: info %d\n", info);
    return 1;
  }

  int i; // C99 compliant.
  for (i = 0; i < nev; ++i) {
    double val = d[i];
    double ref = (N-(nev-1)+i);
    double eps = fabs(val - ref);
    printf("%f - %f - %f\n", val, ref, eps);

    /*eigen value order: smallest -> biggest*/
    if (eps > 1.e-05) return 1;
  }
  return 0;
}

void zMatVec(const a_dcomplex* x, a_dcomplex* y) {
  int i;
  for (i = 0; i < 1000; ++i) y[i] = x[i] * CMPLX(i + 1.0, i + 1.0);
};

int zn() {
  const a_int N      = 1000;
  const a_int nev    = 9;
  const a_int ncv    = 2 * nev + 1;
  const a_int ldv    = N;
  const a_int ldz    = N;
  const a_int lworkl = ncv * (3 * ncv + 5);
  const a_int rvec   = 0;                 // eigenvectors omitted

  const double tol = 0.000001;            // small tol => more stable checks after EV computation.
  const a_dcomplex sigma = CMPLX(0., 0.); // not referenced in this mode

  a_dcomplex resid[N];
  a_dcomplex V[ldv * ncv];
  a_dcomplex z[ldz * nev];
  a_dcomplex d[nev];
  a_dcomplex workd[3 * N];
  a_dcomplex workl[lworkl];
  a_dcomplex workev[2 * ncv];
  double rwork[ncv];
  a_int select[ncv]; // since HOWMNY = 'A', only used as workspace here

  a_int iparam[11], ipntr[14];
  iparam[0] = 1;       // ishift
  iparam[2] = 10 * N;  // on input: maxit; on output: actual iteration
  iparam[3] = 1;       // NB, only 1 allowed
  iparam[6] = 1;       // mode

  char bmat[]   = "I";
  char which[]  = "LM";
  char howmny[] = "A";
  
  a_int info = 0, ido = 0;
  do {
    znaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, rwork, &info);

    zMatVec(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  } while (ido == 1 || ido == -1);

  // check info and number of ev found by arpack
  if (info < 0 || iparam[4] < nev) {
    printf("Error in naupd: iparam[4] %d, nev %d, info %d\n", iparam[4], nev, info);
    return 1;
  }

  zneupd_c(rvec, howmny, select, d, z, ldz, sigma, workev, bmat, N, which, nev,
           tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, rwork, &info);
  if (info < 0) {
    printf("Error in neupd: info %d\n", info);
    return 1;
  }
           
  int i; // C99 compliant.
  for (i = 0; i < nev; ++i) {
    double rval = creal(d[i]);
    double rref = (N-(nev-1)+i);
    double reps = fabs(rval - rref);
    double ival = cimag(d[i]);
    double iref = (N-(nev-1)+i);
    double ieps = fabs(ival - iref);
    printf("%f %f - %f %f - %f %f\n", rval, ival, rref, iref, reps, ieps);

    /*eigen value order: smallest -> biggest*/
    if (reps > 1.e-05 || ieps > 1.e-05) return 1;
  }
  
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
