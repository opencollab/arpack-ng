#include <math.h>
#include <stdio.h>

#include "arpack.h"

/* Nonpositive tolerances must select machine precision on every ICB call. */
static int solve(double tol) {
  enum { n = 10, nev = 2, ncv = 6, lworkl = ncv * (ncv + 8) };
  a_int ido = 0, info = 1;
  a_int iparam[11] = {1, 0, 100, 1, 0, 0, 1}, ipntr[11] = {0};
  a_int select[ncv] = {0};
  double resid[n], v[n * ncv] = {0};
  double workd[3 * n] = {0}, workl[lworkl] = {0};
  double d[nev] = {0}, z[n * nev] = {0};
  int i, j;

  /* A = diag(1,...,10), with a deterministic starting vector. */
  for (i = 0; i < n; ++i) resid[i] = i + 1;
  do {
    dsaupd_c(&ido, "I", n, "LM", nev, tol, resid, ncv, v, n, iparam,
             ipntr, workd, workl, lworkl, &info);
    if (ido == -1 || ido == 1) {
      for (i = 0; i < n; ++i)
        workd[ipntr[1] - 1 + i] = (i + 1) * workd[ipntr[0] - 1 + i];
    } else if (ido != 99) {
      printf("Unexpected reverse communication request: %d\n", (int)ido);
      return 1;
    }
  } while (ido != 99);

  if (info != 0 || iparam[4] != nev) {
    printf("tol=%g: dsaupd info=%d, nconv=%d\n", tol, (int)info,
           (int)iparam[4]);
    return 1;
  }

  dseupd_c(1, "A", select, d, z, n, 0.0, "I", n, "LM", nev, tol,
           resid, ncv, v, n, iparam, ipntr, workd, workl, lworkl, &info);
  if (info != 0) {
    printf("tol=%g: dseupd info=%d\n", tol, (int)info);
    return 1;
  }

  for (j = 0; j < nev; ++j) {
    double residual2 = 0.0, norm2 = 0.0;
    for (i = 0; i < n; ++i) {
      double x = z[n * j + i], r = ((i + 1) - d[j]) * x;
      residual2 += r * r;
      norm2 += x * x;
    }
    if (!(fabs(d[j] - (9 + j)) <= 1.e-10) || !(norm2 >= 0.9) ||
        !(residual2 <= 1.e-20 * norm2)) {
      printf("tol=%g: incorrect eigenpair %d, eigenvalue=%g\n", tol, j, d[j]);
      return 1;
    }
  }
  return 0;
}

int main(void) {
  int result = solve(1.e-12);
  result |= solve(0.0);
  result |= solve(-1.0);
  return result;
}
