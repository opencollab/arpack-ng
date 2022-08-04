/*
 * This example demonstrates the use of ISO_C_BINDING to call parpack
 * (portability). IMPORTANT: MPI communicators MUST be passed from C to Fortran
 * using MPI_Comm_c2f.
 *
 * Just use arpack as you would have normally done but use [ae]upd instead
 * of *[ae]upd_. The main advantage is that checks of the arguments are
 * performed at compile time. Note: to debug parpack, call debug_c.
 * This test program solves for the 9 eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "parpack.hpp"
#include "debug_c.hpp"  // debug parpack.
#include "stat_c.hpp"   // arpack statistics.

template <typename Real>
void diagonal_matrix_vector_product(const Real* x, Real* y) {
  for (int i = 0; i < 1000; ++i) {
    y[i] = static_cast<float>(i + 1) * x[i];
  }
}

template <typename Real>
void real_symmetric_runner(double const& tol_check, arpack::which const& ritz_option) {
  const a_int N      = 1000;
  const a_int nev    = 9;
  const a_int ncv    = 2 * nev + 1;
  const a_int ldv    = N;
  const a_int ldz    = N;
  const a_int lworkl = ncv * (ncv + 8);
  const a_int rvec   = 1;     // need eigenvectors

  const Real tol = 0.000001; // small tol => more stable checks after EV computation.
  const Real sigma = 0.0f;   // not referenced in this mode

  std::vector<Real> resid(N);
  std::vector<Real> V(ldv * ncv);
  std::vector<Real> z(ldz * nev);
  std::vector<Real> d(nev);
  std::vector<Real> workd(3 * N);
  std::vector<Real> workl(lworkl);
  std::vector<a_int> select(ncv); // since HOWMNY = 'A', only used as workspace here

  a_int iparam[11], ipntr[11];
  iparam[0] = 1;      // ishift
  iparam[2] = 10 * N; // on input: maxit; on output: actual iteration
  iparam[3] = 1;      // NB, only 1 allowed
  iparam[6] = 1;      // mode

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);

  a_int info = 0, ido = 0;
  do {
    arpack::saupd(MCW, ido, arpack::bmat::identity, N,
                  ritz_option, nev, tol, resid.data(), ncv,
                  V.data(), ldv, iparam, ipntr, workd.data(),
                  workl.data(), lworkl, info);

    diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  } while  (ido == 1 || ido == -1);

  // check info and number of ev found by arpack.
  if (info < 0 || iparam[4] < nev) { /*arpack may succeed to compute more EV than expected*/
    std::cout << "ERROR in saupd: iparam[4] " << iparam[4] << ", nev " << nev
              << ", info " << info << std::endl;
    throw std::domain_error("Error inside ARPACK routines");
  }

  arpack::seupd(MCW, rvec, arpack::howmny::ritz_vectors, select.data(),
                d.data(), z.data(), ldz, sigma, arpack::bmat::identity, N,
                ritz_option, nev, tol, resid.data(), ncv,
                V.data(), ldv, iparam, ipntr, workd.data(),
                workl.data(), lworkl, info);
  if (info < 0) throw std::runtime_error("Error in seupd, info " + std::to_string(info));

  for (int i = 0; i < nev; ++i) {
    Real val = d[i];
    Real ref = (N - (nev - 1) + i);
    Real eps = std::fabs(val - ref);
    std::cout << "rank " << rank  << " : " << val << " - " << ref << " - " << eps << std::endl;

    /*eigen value order: smallest -> biggest*/
    if (eps > tol_check) throw std::domain_error("Correct eigenvalues not computed");
  }
  std::cout << "------" << std::endl;
}

template <typename Real>
void diagonal_matrix_vector_product(const std::complex<Real>* x, std::complex<Real>* y) {
  for (int i = 0; i < 1000; ++i) {
    // Use complex matrix (i, -i) instead of (i, i): this way "largest_magnitude"
    // and "largest_imaginary" options produce different results that can be checked.
    y[i] = x[i] * std::complex<Real>{Real(i + 1), -Real(i + 1)};
  }
}

template <typename Real>
void complex_nonsymmetric_runner(double const& tol_check, arpack::which const& ritz_option) {
  const a_int N = 1000;
  const a_int nev = 9;
  const a_int ncv = 2 * nev + 1;
  const a_int ldv = N;
  const a_int ldz = N;
  const a_int lworkl = ncv * (3 * ncv + 5);
  const a_int rvec = 0;                       // eigenvectors omitted

  const Real tol = 0.000001;                  // small tol => more stable checks after EV computation.
  const std::complex<Real> sigma(0.0f, 0.0f); // not referenced in this mode

  std::vector<std::complex<Real>> resid(N);
  std::vector<std::complex<Real>> V(ldv * ncv);
  std::vector<std::complex<Real>> z(ldz * nev);
  std::vector<std::complex<Real>> d(nev);
  std::vector<std::complex<Real>> workd(3 * N);
  std::vector<std::complex<Real>> workl(lworkl);
  std::vector<std::complex<Real>> workev(2 * ncv);
  std::vector<Real> rwork(ncv);
  std::vector<a_int> select(ncv); // since HOWMNY = 'A', only used as workspace here

  a_int iparam[11], ipntr[14];
  iparam[0] = 1;       // ishift
  iparam[2] = 10 * N;  // on input: maxit; on output: actual iteration
  iparam[3] = 1;       // NB, only 1 allowed
  iparam[6] = 1;       // mode

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);

  a_int info = 0, ido = 0;
  do {
    arpack::naupd(MCW, ido, arpack::bmat::identity, N,
                  ritz_option, nev, tol, resid.data(), ncv,
                  V.data(), ldv, iparam, ipntr, workd.data(),
                  workl.data(), lworkl, rwork.data(), info);

    diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  } while (ido == 1 || ido == -1);

  // check info and number of ev found by arpack
  if (info < 0 || iparam[4] < nev) { /*arpack may succeed to compute more EV than expected*/
    std::cout << "ERROR in naupd: iparam[4] " << iparam[4] << ", nev " << nev
              << ", info " << info << std::endl;
    throw std::domain_error("Error inside ARPACK routines");
  }

  arpack::neupd(MCW, rvec, arpack::howmny::ritz_vectors, select.data(),
                d.data(), z.data(), ldz, sigma, workev.data(),
                arpack::bmat::identity, N, ritz_option,
                nev, tol, resid.data(), ncv, V.data(), ldv, iparam,
                ipntr, workd.data(), workl.data(), lworkl, rwork.data(), info);
  if (info < 0) throw std::runtime_error("Error in neupd, info " + std::to_string(info));

  if (ritz_option == arpack::which::largest_magnitude) {
    for (int i = 0; i < nev; ++i) {
      Real rval = std::real(d[i]);
      Real rref = static_cast<Real>(N - (nev - 1) + i);
      Real reps = std::fabs(rval - rref);
      Real ival = std::imag(d[i]);
      Real iref = -static_cast<Real>(N - (nev - 1) + i);
      Real ieps = std::fabs(ival - iref);
      std::cout << rval << " " << ival << " - " << rref << " " << iref << " - " << reps << " " << ieps << std::endl;

      if (reps > tol_check || ieps > tol_check) throw std::domain_error("Correct eigenvalues not computed");
    }
  } else if (ritz_option == arpack::which::largest_imaginary) {
    for (int i = 0; i < nev; ++i) {
      Real rval = std::real(d[i]);
      Real rref = static_cast<Real>(nev - i);
      Real reps = std::fabs(rval - rref);
      Real ival = std::imag(d[i]);
      Real iref = -static_cast<Real>(nev - i);
      Real ieps = std::fabs(ival - iref);
      std::cout << rval << " " << ival << " - " << rref << " " << iref << " - " << reps << " " << ieps << std::endl;

      if (reps > tol_check || ieps > tol_check) throw std::domain_error("Correct eigenvalues not computed");
    }
  } else {
    throw std::domain_error("The input Ritz option is not allowed in this test file.");
  }
  std::cout << "------" << std::endl;
}

int main() {
  MPI_Init(NULL, NULL);

  sstats_c();

  try {
    // parpack without debug
    real_symmetric_runner<float>(1., arpack::which::largest_magnitude);
    real_symmetric_runner<float>(1., arpack::which::largest_algebraic);
    real_symmetric_runner<double>(1.e-05, arpack::which::largest_magnitude);
    real_symmetric_runner<double>(1.e-05, arpack::which::largest_algebraic);
  } catch (std::domain_error& e) {
    std::cout << e.what() << std::endl;
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
    complex_nonsymmetric_runner<float>(1., arpack::which::largest_magnitude);
    complex_nonsymmetric_runner<float>(1., arpack::which::largest_imaginary);
    complex_nonsymmetric_runner<double>(1.e-05, arpack::which::largest_magnitude);
    complex_nonsymmetric_runner<double>(1.e-05, arpack::which::largest_imaginary);
  } catch (std::domain_error& e) {
    std::cout << e.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;
}
