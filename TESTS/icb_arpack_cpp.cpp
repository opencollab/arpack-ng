/*
 * This example demonstrates the use of ISO_C_BINDING to call arpack (portability).
 *
 * Just use arpack as you would have normally done, but, use *[ae]upd_c instead of *[ae]upd_.
 * The main advantage is that compiler checks (arguments) are performed at build time.
 * Note: to debug arpack, call debug_c.
 */

#include "arpack.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "debug_c.hpp" // debug arpack.
#include "stat_c.hpp"  // arpack statistics.

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 */
#ifndef BLASINT
#define BLASINT int
#endif

void diagonal_matrix_vector_product(float const* const x, float* const y)
{
    for (int i = 0; i < 1000; ++i)
    {
        y[i] = static_cast<float>(i + 1) * x[i];
    }
}

void real_symmetric_runner()
{
    BLASINT const N = 1000;
    BLASINT const nev = 9;

    BLASINT const ncv = 2 * nev + 1;
    BLASINT const ldv = N;

    BLASINT const ldz = N + 1;

    BLASINT const lworkl = 3 * (ncv * ncv) + 6 * ncv;

    float const tol = 0.0f;
    float const sigma = 0.0f;

    bool const rvec = true;

    std::vector<float> resid(N);
    std::vector<float> V(ncv * N);
    std::vector<float> workd(3 * N, 0.0f);
    std::vector<float> workl(lworkl, 0.0f);
    std::vector<float> d((nev + 1));
    std::vector<float> z((N + 1) * (nev + 1));

    std::array<BLASINT, 11> iparam{};

    iparam[0] = 1;
    iparam[2] = 10 * N;
    iparam[3] = 1;
    iparam[4] = 0; // number of ev found by arpack.
    iparam[6] = 1;

    std::array<BLASINT, 14> ipntr{};

    BLASINT info = 0, ido = 0;

    while (ido != 99)
    {
        /* call arpack like you would have, but, use ssaupd_c instead of ssaupd_ */
        arpack::saupd_c(ido, arpack::bmat::identity, N, arpack::which::largest_magnitude, nev, tol,
                        resid.data(), ncv, V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                        workl.data(), lworkl, info);

        diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
    }

    // check number of ev found by arpack.
    if (iparam[4] != nev || info != 0)
    {
        throw std::domain_error("Error inside ARPACK routines");
    }

    std::vector<int> select(ncv);

    /* call arpack like you would have, but, use sseupd_c instead of sseupd_ */
    arpack::seupd_c(rvec, arpack::howmny::ritz_vectors, select.data(), d.data(), z.data(), ldz,
                    sigma, arpack::bmat::identity, N, arpack::which::largest_magnitude, nev, tol,
                    resid.data(), ncv, V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                    workl.data(), lworkl, info);

    for (int i = 0; i < nev; ++i)
    {
        std::cout << d[i] << "\n";

        if (std::abs(d[i] - static_cast<float>(1000 - (nev - 1) + i)) > 1e-1)
        {
            throw std::domain_error("Correct eigenvalues not computed");
        }
    }
    std::cout << "------\n";
}

void diagonal_matrix_vector_product(std::complex<float> const* const x, std::complex<float>* const y)
{
    for (int i = 0; i < 1000; ++i)
    {
        y[i] = x[i] * std::complex<float>{i + 1.0f, i + 1.0f};
    }
}

void complex_symmetric_runner()
{
    BLASINT const N = 1000;
    BLASINT const nev = 9;

    BLASINT const ncv = 2 * nev + 1;
    BLASINT const ldv = N;

    BLASINT const ldz = N + 1;

    BLASINT const lworkl = 3 * (ncv * ncv) + 6 * ncv;

    float const tol = 0.0f;
    float const sigma = 0.0f;

    bool const rvec = true;

    std::vector<std::complex<float>> resid(N);
    std::vector<std::complex<float>> V(ncv * N);
    std::vector<std::complex<float>> workd(3 * N);
    std::vector<std::complex<float>> workl(lworkl);
    std::vector<std::complex<float>> d(nev + 1);
    std::vector<std::complex<float>> z((N + 1) * (nev + 1));
    std::vector<std::complex<float>> rwork(ncv);
    std::vector<std::complex<float>> workev(2 * ncv);

    std::array<BLASINT, 11> iparam{};
    iparam[0] = 1;
    iparam[2] = 10 * N;
    iparam[3] = 1;
    iparam[4] = 0; // number of ev found by arpack.
    iparam[6] = 1;

    std::array<BLASINT, 14> ipntr{};

    BLASINT info = 0, ido = 0;

    while (ido != 99)
    {
        /* call arpack like you would have, but, use cnaupd_c instead of cnaupd_ */
        arpack::naupd_c(ido, arpack::bmat::identity, N, arpack::which::largest_magnitude, nev, tol,
                        resid.data(), ncv, V.data(), ldv, iparam.data(), ipntr.data(), workd.data(),
                        workl.data(), lworkl, rwork.data(), info);

        diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
    }

    // check number of ev found by arpack.
    if (iparam[4] != nev || info != 0)
    {
        throw std::domain_error("Error inside ARPACK routines");
    }

    std::vector<int> select(ncv);

    /* call arpack like you would have, but, use cneupd_c instead of cneupd_ */
    arpack::neupd_c(rvec, arpack::howmny::ritz_vectors, select.data(), d.data(), z.data(), ldz, sigma,
                    workev.data(), arpack::bmat::identity, N, arpack::which::largest_magnitude, nev,
                    tol, resid.data(), ncv, V.data(), ldv, iparam.data(), ipntr.data(),
                    workd.data(), workl.data(), lworkl, rwork.data(), info);

    for (int i = 0; i < nev; ++i)
    {
        std::cout << d[i] << "\n";

        if (std::abs(std::real(d[i]) - static_cast<float>(1000 - i)) > 1e-1
            || std::abs(std::imag(d[i]) - static_cast<float>(1000 - i)) > 1e-1)
        {
            throw std::domain_error("Correct eigenvalues not computed");
        }
    }
}

int main()
{
    sstats_c();

    real_symmetric_runner(); // arpack without debug.

    int nopx_c, nbx_c, nrorth_c, nitref_c, nrstrt_c;
    float tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c, tsgets_c, tsapps_c, tsconv_c;
    float tnaupd_c, tnaup2_c, tnaitr_c, tneigt_c, tngets_c, tnapps_c, tnconv_c;
    float tcaupd_c, tcaup2_c, tcaitr_c, tceigt_c, tcgets_c, tcapps_c, tcconv_c;
    float tmvopx_c, tmvbx_c, tgetv0_c, titref_c, trvec_c;
    stat_c(nopx_c, nbx_c, nrorth_c, nitref_c, nrstrt_c, tsaupd_c, tsaup2_c, tsaitr_c, tseigt_c,
           tsgets_c, tsapps_c, tsconv_c, tnaupd_c, tnaup2_c, tnaitr_c, tneigt_c, tngets_c, tnapps_c,
           tnconv_c, tcaupd_c, tcaup2_c, tcaitr_c, tceigt_c, tcgets_c, tcapps_c, tcconv_c, tmvopx_c,
           tmvbx_c, tgetv0_c, titref_c, trvec_c);
    std::cout << "Timers : nopx " << nopx_c << ", tmvopx " << tmvopx_c;
    std::cout << " - nbx " << nbx_c << ", tmvbx " << tmvbx_c << std::endl;

    std::cout << "------" << std::endl;

    debug_c(6, -6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1); // set debug flags.

    complex_symmetric_runner(); // arpack with debug.

    return 0;
}
