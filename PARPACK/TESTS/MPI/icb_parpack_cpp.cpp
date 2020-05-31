/*
 * This example demonstrates the use of ISO_C_BINDING to call parpack
 * (portability). IMPORTANT: MPI communicators MUST be passed from C to Fortran
 * using MPI_Comm_c2f.
 *
 * Just use arpack as you would have normally done but use [ae]upd instead
 * of *[ae]upd_. The main advantage is that checks of the arguments are
 * performed at compile time. Note: to debug parpack, call debug_c.
 * This test program solves for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "debug_c.hpp"  // debug parpack.
#include "parpack.hpp"
#include "stat_c.hpp"  // arpack statistics.

void diagonal_matrix_vector_product(float const* const x, float* const y) {
  for (int i = 0; i < 1000; ++i) {
    y[i] = static_cast<float>(i + 1) * x[i];
  }
}

void real_symmetric_runner() {
  a_int N = 1000;
  a_int nev = 3;
  a_int ncv = 2 * nev + 1;
  a_int ldz = N + 1;
  a_int lworkl = 3 * (ncv * ncv) + 6 * ncv;
  a_int ldv = N;

  a_int rvec = 1;
  float tol = 0.000001; // small tol => more stable checks after EV computation.
  float sigma = 0.0f;

  std::array<a_int, 14> ipntr;

  std::vector<float> workd(3 * N, 0.0f);
  std::vector<float> workl(3 * (ncv * ncv) + 6 * ncv, 0.0f);
  std::vector<float> V(ncv * N);
  std::vector<float> d(nev + 1);
  std::vector<float> z((N + 1) * (nev + 1));
  std::vector<float> resid(N);
  std::vector<a_int> select(ncv);
  for (int i = 0; i < ncv; i++) select[i] = 1;

  a_int info = 0;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<a_int, 11> iparam;
  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);

  a_int ido = 0;

  while (ido != 99) {
    arpack::saupd(MCW, ido, arpack::bmat::identity, N,
                  arpack::which::largest_magnitude, nev, tol, resid.data(), ncv,
                  V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                  workl.data(), lworkl, info);

    diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]),
                                   &(workd[ipntr[1] - 1]));
  }
  // check number of ev found by arpack.
  if (iparam[4] < nev /*arpack may succeed to compute more EV than expected*/ ||
      info != 0) {
    std::cout << "ERROR: iparam[4] " << iparam[4] << ", nev " << nev
              << ", info " << info << std::endl;
    throw std::domain_error("Error inside ARPACK routines");
  }

  arpack::seupd(MCW, rvec, arpack::howmny::ritz_vectors, select.data(),
                d.data(), z.data(), ldz, sigma, arpack::bmat::identity, N,
                arpack::which::largest_magnitude, nev, tol, resid.data(), ncv,
                V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                workl.data(), lworkl, info);

  for (int i = 0; i < nev; ++i) {
    float val = d[i];
    float ref = (N - (nev - 1) + i);
    float eps = std::fabs(val - ref);
    std::cout << "rank " << rank  << " : " << val << " - " << ref << " - " << eps << std::endl;

    /*eigen value order: smallest -> biggest*/
    if (eps > 1.) {
      throw std::domain_error("Correct eigenvalues not computed");
    }
  }
}

void diagonal_matrix_vector_product(std::complex<float> const* const x,
                                    std::complex<float>* const y) {
  for (int i = 0; i < 1000; ++i) {
    y[i] = x[i] * std::complex<float>{i + 1.0f, i + 1.0f};
  }
}

void complex_symmetric_runner() {
  a_int N = 1000;
  a_int nev = 1;
  a_int ncv = 2 * nev + 1;
  a_int ldv = N;
  a_int ldz = N + 1;

  float tol = 0.000001; // small tol => more stable checks after EV computation.
  a_int rvec = 0;
  std::complex<float> sigma(0.0f, 0.0f);

  std::vector<std::complex<float>> resid(N);
  std::vector<std::complex<float>> V(ncv * N);
  std::vector<std::complex<float>> workd(3 * N);
  std::vector<std::complex<float>> d(nev + 1);
  std::vector<std::complex<float>> z((N + 1) * (nev + 1));
  std::vector<a_int> select(ncv);
  for (int i = 0; i < ncv; i++) select[i] = 1;

  a_int lworkl = 3 * (ncv * ncv) + 6 * ncv;
  std::vector<std::complex<float>> workl(lworkl);

  std::vector<std::complex<float>> rwork(ncv);
  std::vector<std::complex<float>> workev(2 * ncv);

  a_int info = 0;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<a_int, 11> iparam;
  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  std::array<a_int, 14> ipntr;

  MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);

  a_int ido = 0;

  while (ido != 99) {
    arpack::naupd(MCW, ido, arpack::bmat::identity, N,
                  arpack::which::largest_magnitude, nev, tol, resid.data(), ncv,
                  V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                  workl.data(), lworkl, rwork.data(), info);

    diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]),
                                   &(workd[ipntr[1] - 1]));
  }

  // check number of ev found by arpack
  if (iparam[4] < nev /*arpack may succeed to compute more EV than expected*/ ||
      info != 0) {
    std::cout << "ERROR: iparam[4] " << iparam[4] << ", nev " << nev
              << ", info " << info << std::endl;
    throw std::domain_error("Error inside ARPACK routines");
  }

  arpack::neupd(MCW, rvec, arpack::howmny::ritz_vectors, select.data(),
                d.data(), z.data(), ldz, sigma, workev.data(),
                arpack::bmat::identity, N, arpack::which::largest_magnitude,
                nev, tol, resid.data(), ncv, V.data(), ldv, iparam.data(),
                ipntr.data(), workd.data(), workl.data(), lworkl, rwork.data(),
                info);

  for (int i = 0; i < nev; ++i) {
    float rval = std::real(d[i]);
    float rref = (N-(nev-1)+i);
    float reps = std::fabs(rval - rref);
    float ival = std::imag(d[i]);
    float iref = (N-(nev-1)+i);
    float ieps = std::fabs(ival - iref);
    std::cout << "rank " << rank << " : " << rval << " " << ival << " - " << rref << " " << iref << " - " << reps << " " << ieps << std::endl;

    /*eigen value order: smallest -> biggest*/
    if (reps > 1. || ieps > 1.) {
      throw std::domain_error("Correct eigenvalues not computed");
    }
  }
}

int main() {
  MPI_Init(NULL, NULL);

  sstats_c();

  try {
    // parpack without debug
    real_symmetric_runner();
  } catch (std::domain_error& e) {
    std::cout << e.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  a_int nopx_c, nbx_c, nrorth_c, nitref_c, nrstrt_c;
  float tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c;
  float tnaupd_c, tnaup2_c, tnaitr_c, tneigt_c, tngets_c, tnapps_c, tnconv_c;
  float tcaupd_c, tcaup2_c, tcaitr_c, tceigt_c, tcgets_c, tcapps_c, tcconv_c;
  float tmvopx_c, tmvbx_c, tgetv0_c, titref_c, trvec_c;

  stat_c(nopx_c, nbx_c, nrorth_c, nitref_c, nrstrt_c, tsaupd_c, tsaup2_c,
         tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c, tnaupd_c, tnaup2_c,
         tnaitr_c, tneigt_c, tngets_c, tnapps_c, tnconv_c, tcaupd_c, tcaup2_c,
         tcaitr_c, tceigt_c, tcgets_c, tcapps_c, tcconv_c, tmvopx_c, tmvbx_c,
         tgetv0_c, titref_c, trvec_c);

  std::cout << "Timers : nopx " << nopx_c << ", tmvopx " << tmvopx_c;
  std::cout << " - nbx " << nbx_c << ", tmvbx " << tmvbx_c << std::endl;

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) std::cout << "------" << std::endl;

  // set debug flags.
  debug_c(6, -6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1);

  try {
    complex_symmetric_runner();
  } catch (std::domain_error& e) {
    std::cout << e.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;
}
