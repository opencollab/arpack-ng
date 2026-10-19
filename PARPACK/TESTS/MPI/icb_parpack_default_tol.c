#include <math.h>
#include <stdio.h>

#include "parpack.h"

/* Nonpositive tolerances must select machine precision on every ICB call. */
enum { n = 10, nev = 2, ncv = 6 };

static MPI_Fint comm;
static int rank, nprocs, start, nloc;

/* y = A*x on the local rows: A has diagonal 1,...,n and, if upper is set,
 * superdiagonal 1. */
static void matvec(int upper, const double* x, double* y) {
  double next = 0.0;
  int i;

  if (upper) /* The last local row needs the first entry of the next rank. */
    MPI_Sendrecv(x, 1, MPI_DOUBLE, rank > 0 ? rank - 1 : MPI_PROC_NULL, 0,
                 &next, 1, MPI_DOUBLE,
                 rank + 1 < nprocs ? rank + 1 : MPI_PROC_NULL, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  for (i = 0; i < nloc; ++i) {
    y[i] = (start + i + 1) * x[i];
    if (upper) y[i] += i + 1 < nloc ? x[i + 1] : next;
  }
}

/* Check that z is a unit eigenvector of A with eigenvalue lambda. */
static int check(int upper, double lambda, const double* z) {
  double y[n], sums[2] = {0.0, 0.0};
  int i;

  matvec(upper, z, y);
  for (i = 0; i < nloc; ++i) {
    double r = y[i] - lambda * z[i];
    sums[0] += r * r;
    sums[1] += z[i] * z[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return !(sums[1] >= 0.9) || !(sums[0] <= 1.e-20 * sums[1]);
}

static int solve(double tol) {
  enum { lworkl = ncv * (ncv + 8) };
  a_int ido = 0, info = 1;
  a_int iparam[11] = {1, 0, 100, 1, 0, 0, 1}, ipntr[11] = {0};
  a_int select[ncv] = {0};
  double resid[n], v[n * ncv] = {0};
  double workd[3 * n] = {0}, workl[lworkl] = {0};
  double d[nev] = {0}, z[n * nev] = {0};
  int i, j;

  /* A = diag(1,...,n), with a deterministic starting vector. */
  for (i = 0; i < nloc; ++i) resid[i] = start + i + 1;
  do {
    pdsaupd_c(comm, &ido, "I", nloc, "LM", nev, tol, resid, ncv, v, nloc,
              iparam, ipntr, workd, workl, lworkl, &info);
    if (ido == -1 || ido == 1) {
      matvec(0, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
    } else if (ido != 99) {
      printf("Unexpected reverse communication request: %d\n", (int)ido);
      return 1;
    }
  } while (ido != 99);

  if (info != 0 || iparam[4] != nev) {
    printf("tol=%g: pdsaupd info=%d, nconv=%d\n", tol, (int)info,
           (int)iparam[4]);
    return 1;
  }

  pdseupd_c(comm, 1, "A", select, d, z, nloc, 0.0, "I", nloc, "LM", nev, tol,
            resid, ncv, v, nloc, iparam, ipntr, workd, workl, lworkl, &info);
  if (info != 0) {
    printf("tol=%g: pdseupd info=%d\n", tol, (int)info);
    return 1;
  }

  for (j = 0; j < nev; ++j) {
    if (check(0, d[j], &z[nloc * j]) | !(fabs(d[j] - (n - 1 + j)) <= 1.e-10)) {
      printf("tol=%g: incorrect eigenpair %d, eigenvalue=%g\n", tol, j, d[j]);
      return 1;
    }
  }
  return 0;
}

static int solve_nonsymmetric(double tol) {
  enum { lworkl = 3 * ncv * (ncv + 2) };
  a_int ido = 0, info = 1;
  a_int iparam[11] = {1, 0, 100, 1, 0, 0, 1}, ipntr[14] = {0};
  a_int select[ncv] = {0};
  double resid[n], v[n * ncv] = {0};
  double workd[3 * n] = {0}, workl[lworkl] = {0}, workev[3 * ncv] = {0};
  double dr[nev + 1] = {0}, di[nev + 1] = {0}, z[n * (nev + 1)] = {0};
  int i, j;

  /* Upper bidiagonal A has diagonal 1,...,n and superdiagonal 1. */
  for (i = 0; i < nloc; ++i) resid[i] = start + i + 1;
  do {
    pdnaupd_c(comm, &ido, "I", nloc, "LM", nev, tol, resid, ncv, v, nloc,
              iparam, ipntr, workd, workl, lworkl, &info);
    if (ido == -1 || ido == 1) {
      matvec(1, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
    } else if (ido != 99) {
      printf("Unexpected reverse communication request: %d\n", (int)ido);
      return 1;
    }
  } while (ido != 99);

  if (info != 0 || iparam[4] != nev) {
    printf("tol=%g: pdnaupd info=%d, nconv=%d\n", tol, (int)info,
           (int)iparam[4]);
    return 1;
  }

  pdneupd_c(comm, 1, "A", select, dr, di, z, nloc, 0.0, 0.0, workev, "I",
            nloc, "LM", nev, tol, resid, ncv, v, nloc, iparam, ipntr, workd,
            workl, lworkl, &info);
  if (info != 0) {
    printf("tol=%g: pdneupd info=%d\n", tol, (int)info);
    return 1;
  }

  for (j = 0; j < nev; ++j) {
    double expected = dr[0] < dr[1] ? n - 1 + j : n - j;
    if (check(1, dr[j], &z[nloc * j]) | !(fabs(dr[j] - expected) <= 1.e-10) ||
        !(fabs(di[j]) <= 1.e-10)) {
      printf("tol=%g: incorrect nonsymmetric eigenpair %d, eigenvalue=%g+%gi\n",
             tol, j, dr[j], di[j]);
      return 1;
    }
  }
  return 0;
}

int main(void) {
  int result;

  MPI_Init(NULL, NULL);
  comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  nloc = n / nprocs + (rank < n % nprocs);
  start = rank * (n / nprocs) + (rank < n % nprocs ? rank : n % nprocs);

  result = solve(1.e-12);
  result |= solve(0.0);
  result |= solve(-1.0);
  result |= solve_nonsymmetric(0.0);
  result |= solve_nonsymmetric(-1.0);

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return result;
}
