#ifndef __ARPACK_HPP__
#define __ARPACK_HPP__

#include <complex.h>
#include <complex>

namespace arpack {
enum class which : int {
  /// 'LA' - compute the NEV largest (algebraic) eigenvalues
  largest_algebraic,
  /// 'SA' - compute the NEV smallest (algebraic) eigenvalues.
  smallest_algebraic,
  /// 'LM' - compute the NEV largest (in magnitude) eigenvalues.
  largest_magnitude,
  /// 'SM' - compute the NEV smallest (in magnitude) eigenvalues.
  smallest_magnitude,
  /// 'BE' - compute NEV eigenvalues, half from each end of the
  /// spectrum.  When NEV is odd, compute one more from the
  /// high end than from the low end.
  both_ends
};

enum class bmat : int {
  /// B = 'I' -> standard eigenvalue problem A*x = lambda*x
  identity,
  /// B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
  generalized
};

enum class howmny : int {
  /// 'A' Compute NEV Ritz vectors
  ritz_vectors,
  /// 'P' Compute NEV Schur vectors;
  schur_vectors,
  /// 'S' compute some of the Ritz vectors, specified by the logical array
  /// SELECT.
  ritz_specified
};

namespace internal {
/*
 * From C++, arpack does not exist.
 * Arpack is Fortran. ISO_C_BINDING is a gateway from Fortran to C, not C++.
 * However C++ can interface to C which is why you find C types in the
 * arpack.hpp
 */
extern "C" {
void ssaupd_c(int& ido, const char* bmat, int n, const char* which, int nev,
              float tol, float* resid, int ncv, float* v, int ldv, int* iparam,
              int* ipntr, float* workd, float* workl, int lworkl, int& info);

void sseupd_c(bool rvec, const char* howmny, int* select, float* d, float* z,
              int ldz, float sigma, const char* bmat, int n, const char* which,
              int nev, float tol, float* resid, int ncv, float* v, int ldv,
              int* iparam, int* ipntr, float* workd, float* workl, int lworkl,
              int& info);

void dsaupd_c(int& ido, const char* bmat, int n, const char* which, int nev,
              double tol, double* resid, int ncv, double* v, int ldv,
              int* iparam, int* ipntr, double* workd, double* workl, int lworkl,
              int& info);

void dseupd_c(bool rvec, const char* howmny, int* select, double* d, double* z,
              int ldz, double sigma, const char* bmat, int n, const char* which,
              int nev, double tol, double* resid, int ncv, double* v, int ldv,
              int* iparam, int* ipntr, double* workd, double* workl, int lworkl,
              int& info);

void snaupd_c(int& ido, const char* bmat, int n, const char* which, int nev,
              float tol, float* resid, int ncv, float* v, int ldv, int* iparam,
              int* ipntr, float* workd, float* workl, int lworkl, int& info);

void sneupd_c(bool rvec, const char* howmny, int* select, float* dr, float* di,
              float* z, int ldz, float sigmar, float sigmai, const char* bmat,
              int n, const char* which, int nev, float tol, float* resid,
              int ncv, float* v, int ldv, int* iparam, int* ipntr, float* workd,
              float* workl, int lworkl, int& info);

void dnaupd_c(int& ido, const char* bmat, int n, const char* which, int nev,
              double tol, double* resid, int ncv, double* v, int ldv,
              int* iparam, int* ipntr, double* workd, double* workl, int lworkl,
              int& info);

void dneupd_c(bool rvec, const char* howmny, int* select, double* dr,
              double* di, double* z, int ldz, double sigmar, double sigmai,
              const char* bmat, int n, const char* which, int nev, double tol,
              double* resid, int ncv, double* v, int ldv, int* iparam,
              int* ipntr, double* workd, double* workl, int lworkl, int& info);

void cnaupd_c(int& ido, const char* bmat, int n, const char* which, int nev,
              float tol, float _Complex* resid, int ncv, float _Complex* v,
              int ldv, int* iparam, int* ipntr, float _Complex* workd,
              float _Complex* workl, int lworkl, float _Complex* rwork,
              int& info);

void cneupd_c(bool rvec, const char* howmny, int* select, float _Complex* d,
              float _Complex* z, int ldz, float _Complex const* sigma,
              float _Complex* workev, const char* bmat, int n,
              const char* which, int nev, float tol, float _Complex* resid,
              int ncv, float _Complex* v, int ldv, int* iparam, int* ipntr,
              float _Complex* workd, float _Complex* workl, int lworkl,
              float _Complex* rwork, int& info);

void znaupd_c(int& ido, const char* bmat, int n, const char* which, int nev,
              double tol, double _Complex* resid, int ncv, double _Complex* v,
              int ldv, int* iparam, int* ipntr, double _Complex* workd,
              double _Complex* workl, int lworkl, double _Complex* rwork,
              int& info);

void zneupd_c(bool rvec, const char* howmny, int* select, double _Complex* d,
              double _Complex* z, int ldz, double _Complex const* sigma,
              double _Complex* workev, const char* bmat, int n,
              const char* which, int nev, double tol, double _Complex* resid,
              int ncv, double _Complex* v, int ldv, int* iparam, int* ipntr,
              double _Complex* workd, double _Complex* workl, int lworkl,
              double _Complex* rwork, int& info);
}

inline char const* convert_to_char(which const option) {
  switch (option) {
    case which::largest_algebraic: {
      return "LA";
      break;
    }
    case which::smallest_algebraic: {
      return "SA";
      break;
    }
    case which::largest_magnitude: {
      return "LM";
      break;
    }
    case which::smallest_magnitude: {
      return "SM";
      break;
    }
    case which::both_ends: {
      return "BE";
      break;
    }
  }
  return "LM";
}

inline char const* convert_to_char(bmat const option) {
  return option == bmat::identity ? "I" : "B";
}

inline char const* convert_to_char(howmny const option) {
  switch (option) {
    case howmny::ritz_vectors: {
      return "A";
      break;
    }
    case howmny::schur_vectors: {
      return "P";
      break;
    }
    case howmny::ritz_specified: {
      return "S";
      break;
    }
  }
  return "A";
}
}  // namespace internal

inline void saupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::ssaupd_c(ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void seupd(bool rvec, howmny const howmny_option, int* select, float* d,
                  float* z, int ldz, float sigma, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::sseupd_c(rvec, internal::convert_to_char(howmny_option), select, d,
                     z, ldz, sigma, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void saupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::dsaupd_c(ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void seupd(bool rvec, howmny const howmny_option, int* select, double* d,
                  double* z, int ldz, double sigma, bmat const bmat_option,
                  int n, which const ritz_option, int nev, double tol,
                  double* resid, int ncv, double* v, int ldv, int* iparam,
                  int* ipntr, double* workd, double* workl, int lworkl,
                  int& info) {
  internal::dseupd_c(rvec, internal::convert_to_char(howmny_option), select, d,
                     z, ldz, sigma, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::snaupd_c(ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void neupd(bool rvec, howmny const howmny_option, int* select, float* dr,
                  float* di, float* z, int ldz, float sigmar, float sigmai,
                  bmat const bmat_option, int n, which const ritz_option,
                  int nev, float tol, float* resid, int ncv, float* v, int ldv,
                  int* iparam, int* ipntr, float* workd, float* workl,
                  int lworkl, int& info) {
  internal::sneupd_c(rvec, internal::convert_to_char(howmny_option), select, dr,
                     di, z, ldz, sigmar, sigmai,
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::dnaupd_c(ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void neupd(bool rvec, howmny const howmny_option, int* select,
                  double* dr, double* di, double* z, int ldz, double sigmar,
                  double sigmai, bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::dneupd_c(rvec, internal::convert_to_char(howmny_option), select, dr,
                     di, z, ldz, sigmar, sigmai,
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol, resid,
                     ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, float tol,
                  std::complex<float>* resid, int ncv, std::complex<float>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<float>* workd,
                  std::complex<float>* workl, int lworkl,
                  std::complex<float>* rwork, int& info) {
  internal::cnaupd_c(ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<_Complex float*>(resid), ncv,
                     reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<_Complex float*>(workd),
                     reinterpret_cast<_Complex float*>(workl), lworkl,
                     reinterpret_cast<_Complex float*>(rwork), info);
}

inline void neupd(bool rvec, howmny const howmny_option, int* select,
                  std::complex<float>* d, std::complex<float>* z, int ldz,
                  std::complex<float> sigma, std::complex<float>* workev,
                  bmat const bmat_option, int n, which const ritz_option,
                  int nev, float tol, std::complex<float>* resid, int ncv,
                  std::complex<float>* v, int ldv, int* iparam, int* ipntr,
                  std::complex<float>* workd, std::complex<float>* workl,
                  int lworkl, std::complex<float>* rwork, int& info) {
  const _Complex float csigma=std::real(sigma) + std::imag(sigma) * I;
  internal::cneupd_c(rvec, internal::convert_to_char(howmny_option), select,
                     reinterpret_cast<_Complex float*>(d),
                     reinterpret_cast<_Complex float*>(z), ldz,
                     &csigma,
                     reinterpret_cast<_Complex float*>(workev),
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<_Complex float*>(resid), ncv,
                     reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<_Complex float*>(workd),
                     reinterpret_cast<_Complex float*>(workl), lworkl,
                     reinterpret_cast<_Complex float*>(rwork), info);
}

inline void naupd(int& ido, bmat const bmat_option, int n,
                  which const ritz_option, int nev, double tol,
                  std::complex<double>* resid, int ncv, std::complex<double>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<double>* workd,
                  std::complex<double>* workl, int lworkl,
                  std::complex<double>* rwork, int& info) {
  internal::znaupd_c(ido, internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<_Complex double*>(resid), ncv,
                     reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<_Complex double*>(workd),
                     reinterpret_cast<_Complex double*>(workl), lworkl,
                     reinterpret_cast<_Complex double*>(rwork), info);
}

inline void neupd(bool rvec, howmny const howmny_option, int* select,
                  std::complex<double>* d, std::complex<double>* z, int ldz,
                  std::complex<double> sigma, std::complex<double>* workev,
                  bmat const bmat_option, int n, which const ritz_option,
                  int nev, double tol, std::complex<double>* resid, int ncv,
                  std::complex<double>* v, int ldv, int* iparam, int* ipntr,
                  std::complex<double>* workd, std::complex<double>* workl,
                  int lworkl, std::complex<double>* rwork, int& info) {
  const _Complex double csigma=std::real(sigma) + _Complex_I * std::imag(sigma);
  internal::zneupd_c(rvec, internal::convert_to_char(howmny_option), select,
                     reinterpret_cast<_Complex double*>(d),
                     reinterpret_cast<_Complex double*>(z), ldz,
                     &csigma,
                     reinterpret_cast<_Complex double*>(workev),
                     internal::convert_to_char(bmat_option), n,
                     internal::convert_to_char(ritz_option), nev, tol,
                     reinterpret_cast<_Complex double*>(resid), ncv,
                     reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
                     reinterpret_cast<_Complex double*>(workd),
                     reinterpret_cast<_Complex double*>(workl), lworkl,
                     reinterpret_cast<_Complex double*>(rwork), info);
}
}  // namespace arpack

#endif
