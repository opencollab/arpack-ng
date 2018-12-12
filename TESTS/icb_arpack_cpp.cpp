/*
 * This example demonstrates the use of C++ bindings to call arpack.
 *
 * Use arpack as you would have normally done, but, use [ae]upd instead
 * of *[ae]upd_. The main advantage is that compiler checks the argument types
 * and the correct function is called based on the type (float vs double vs
 * complex). Note: to debug arpack, call debug_c. This is a test program to
 * solve for the 9 largest eigenvalues of A*x = lambda*x where A is the diagonal
 * matrix with entries 1000, 999, ... , 2, 1 on the diagonal.
 */

#include "arpack.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "debug_c.hpp"  // debug arpack.
#include "stat_c.hpp"   // arpack statistics.

template<typename Real>
void diagonal_matrix_vector_product(Real const* const x, Real* const y) {
  for (int i = 0; i < 1000; ++i) {
    y[i] = static_cast<Real>(i + 1) * x[i];
  }
}

template<typename Real>
void real_symmetric_runner() {
  a_int const N = 1000;
  a_int const nev = 9;

  a_int const ncv = 2 * nev + 1;
  a_int const ldv = N;

  a_int const ldz = N + 1;

  a_int const lworkl = 3 * (ncv * ncv) + 6 * ncv;

  Real const tol = 0.0;
  Real const sigma = 0.0;

  bool const rvec = true;

  std::vector<Real> resid(N);
  std::vector<Real> V(ncv * N);
  std::vector<Real> workd(3 * N, 0.0);
  std::vector<Real> workl(lworkl, 0.0);
  std::vector<Real> d((nev + 1));
  std::vector<Real> z((N + 1) * (nev + 1));

  std::array<a_int, 11> iparam{};

  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  std::array<a_int, 14> ipntr{};

  a_int info = 0, ido = 0;

  while (ido != 99) {
    arpack::saupd(ido, arpack::bmat::identity, N,
                  arpack::which::largest_magnitude, nev, tol, resid.data(), ncv,
                  V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                  workl.data(), lworkl, info);

    diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]),
                                   &(workd[ipntr[1] - 1]));
  }

  // check number of ev found by arpack.
  if (iparam[4] != nev || info != 0) {
    throw std::domain_error("Error inside ARPACK routines");
  }

  std::vector<a_int> select(ncv);

  arpack::seupd(rvec, arpack::howmny::ritz_vectors, select.data(), d.data(),
                z.data(), ldz, sigma, arpack::bmat::identity, N,
                arpack::which::largest_magnitude, nev, tol, resid.data(), ncv,
                V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                workl.data(), lworkl, info);

  for (int i = 0; i < nev; ++i) {
    std::cout << d[i] << "\n";

    if (std::abs(d[i] - static_cast<Real>(1000 - (nev - 1) + i)) > 1e-1) {
      throw std::domain_error("Correct eigenvalues not computed");
    }
  }
  std::cout << "------\n";
}

template<typename Real>
void diagonal_matrix_vector_product(std::complex<Real> const* const x,
                                    std::complex<Real>* const y) {
  for (int i = 0; i < 1000; ++i) {
    y[i] = x[i] * std::complex<Real>{Real(i + 1), Real(i + 1)};
  }
}

template<typename Real>
void complex_symmetric_runner() {
  a_int const N = 1000;
  a_int const nev = 9;

  a_int const ncv = 2 * nev + 1;
  a_int const ldv = N;

  a_int const ldz = N + 1;

  a_int const lworkl = 3 * (ncv * ncv) + 6 * ncv;

  Real const tol = 0.0;
  std::complex<Real> const sigma(0.0, 0.0);

  bool const rvec = true;

  std::vector<std::complex<Real>> resid(N);
  std::vector<std::complex<Real>> V(ncv * N);
  std::vector<std::complex<Real>> workd(3 * N);
  std::vector<std::complex<Real>> workl(lworkl);
  std::vector<std::complex<Real>> d(nev + 1);
  std::vector<std::complex<Real>> z((N + 1) * (nev + 1));
  std::vector<Real> rwork(ncv);
  std::vector<std::complex<Real>> workev(2 * ncv);

  std::array<a_int, 11> iparam{};
  iparam[0] = 1;
  iparam[2] = 10 * N;
  iparam[3] = 1;
  iparam[4] = 0;  // number of ev found by arpack.
  iparam[6] = 1;

  std::array<a_int, 14> ipntr{};

  a_int info = 0, ido = 0;

  while (ido != 99) {
    arpack::naupd(ido, arpack::bmat::identity, N,
                  arpack::which::largest_magnitude, nev, tol, resid.data(), ncv,
                  V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                  workl.data(), lworkl, rwork.data(), info);

    diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]),
                                   &(workd[ipntr[1] - 1]));
  }

  // check number of ev found by arpack.
  if (iparam[4] != nev || info != 0) {
    throw std::domain_error("Error inside ARPACK routines");
  }

  std::vector<a_int> select(ncv);

  arpack::neupd(rvec, arpack::howmny::ritz_vectors, select.data(), d.data(),
                z.data(), ldz, sigma, workev.data(), arpack::bmat::identity, N,
                arpack::which::largest_magnitude, nev, tol, resid.data(), ncv,
                V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                workl.data(), lworkl, rwork.data(), info);

  for (int i = 0; i < nev; ++i) {
    std::cout << d[i] << "\n";

    if (std::abs(std::real(d[i]) - static_cast<Real>(1000 - i)) > 1e-1 ||
        std::abs(std::imag(d[i]) - static_cast<Real>(1000 - i)) > 1e-1) {
      throw std::domain_error("Correct eigenvalues not computed");
    }
  }
}

int main() {
  sstats_c();

  // arpack without debug
  real_symmetric_runner<float>();
  real_symmetric_runner<double>();

  int nopx_c, nbx_c, nrorth_c, nitref_c, nrstrt_c;
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

  std::cout << "------" << std::endl;

  // set debug flags
  debug_c(6, -6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1);

  // arpack with debug
  complex_symmetric_runner<float>();
  complex_symmetric_runner<double>();

  return 0;
}
