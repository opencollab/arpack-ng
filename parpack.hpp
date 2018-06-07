
#ifndef __PARPACK_HPP__
#define __PARPACK_HPP__

#include "arpack.hpp"

#include <mpi.h>

#include <complex.h>
#include <complex>

namespace arpack {

namespace internal {
#include "parpack.h"
}  // namespace internal

inline void saupd(MPI_Fint comm, int& ido, bmat const bmat_option, int n,
                  which const which_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::pssaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void seupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  int* select, float* d, float* z, int ldz, float sigma,
                  bmat const bmat_option, int n, which const which_option,
                  int nev, float tol, float* resid, int ncv, float* v, int ldv,
                  int* iparam, int* ipntr, float* workd, float* workl,
                  int lworkl, int& info) {
  internal::psseupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, d, z, ldz, sigma,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void saupd(MPI_Fint comm, int& ido, bmat const bmat_option, int n,
                  which const which_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::pdsaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void seupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  int* select, double* d, double* z, int ldz, double sigma,
                  bmat const bmat_option, int n, which const which_option,
                  int nev, double tol, double* resid, int ncv, double* v,
                  int ldv, int* iparam, int* ipntr, double* workd,
                  double* workl, int lworkl, int& info) {
  internal::pdseupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, d, z, ldz, sigma,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(MPI_Fint comm, int& ido, bmat const bmat_option, int n,
                  which const which_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::psnaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  int* select, float* dr, float* di, float* z, int ldz,
                  float sigmar, float sigmai, bmat const bmat_option, int n,
                  which const which_option, int nev, float tol, float* resid,
                  int ncv, float* v, int ldv, int* iparam, int* ipntr,
                  float* workd, float* workl, int lworkl, int& info) {
  internal::psneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, dr, di, z, ldz, sigmar, sigmai,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(MPI_Fint comm, int& ido, bmat const bmat_option, int n,
                  which const which_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::pdnaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  int* select, double* dr, double* di, double* z, int ldz,
                  double sigmar, double sigmai, bmat const bmat_option, int n,
                  which const which_option, int nev, double tol, double* resid,
                  int ncv, double* v, int ldv, int* iparam, int* ipntr,
                  double* workd, double* workl, int lworkl, int& info) {
  internal::pdneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, dr, di, z, ldz, sigmar, sigmai,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(MPI_Fint comm, int& ido, bmat const bmat_option, int n,
                  which const which_option, int nev, float tol,
                  std::complex<float>* resid, int ncv, std::complex<float>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<float>* workd,
                  std::complex<float>* workl, int lworkl,
                  std::complex<float>* rwork, int& info) {
  internal::pcnaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex float*>(resid), ncv,
                      reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex float*>(workd),
                      reinterpret_cast<_Complex float*>(workl), lworkl,
                      reinterpret_cast<_Complex float*>(rwork), &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  int* select, std::complex<float>* d, std::complex<float>* z,
                  int ldz, std::complex<float> sigma,
                  std::complex<float>* workev, bmat const bmat_option, int n,
                  which const which_option, int nev, float tol,
                  std::complex<float>* resid, int ncv, std::complex<float>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<float>* workd,
                  std::complex<float>* workl, int lworkl,
                  std::complex<float>* rwork, int& info)

{
  internal::pcneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, reinterpret_cast<_Complex float*>(d),
                      reinterpret_cast<_Complex float*>(z), ldz,
                      std::real(sigma) + _Complex_I * std::imag(sigma),
                      reinterpret_cast<_Complex float*>(workev),
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex float*>(resid), ncv,
                      reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex float*>(workd),
                      reinterpret_cast<_Complex float*>(workl), lworkl,
                      reinterpret_cast<_Complex float*>(rwork), &info);
}

inline void naupd(MPI_Fint comm, int& ido, bmat const bmat_option, int n,
                  which const which_option, int nev, double tol,
                  std::complex<double>* resid, int ncv, std::complex<double>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<double>* workd,
                  std::complex<double>* workl, int lworkl,
                  std::complex<double>* rwork, int& info) {
  internal::pznaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex double*>(resid), ncv,
                      reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex double*>(workd),
                      reinterpret_cast<_Complex double*>(workl), lworkl,
                      reinterpret_cast<_Complex double*>(rwork), &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  int* select, std::complex<double>* d, std::complex<double>* z,
                  int ldz, std::complex<double> sigma,
                  std::complex<double>* workev, bmat const bmat_option, int n,
                  which const which_option, int nev, double tol,
                  std::complex<double>* resid, int ncv, std::complex<double>* v,
                  int ldv, int* iparam, int* ipntr, std::complex<double>* workd,
                  std::complex<double>* workl, int lworkl,
                  std::complex<double>* rwork, int& info) {
  internal::pzneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, reinterpret_cast<_Complex double*>(d),
                      reinterpret_cast<_Complex double*>(z), ldz,
                      std::real(sigma) + _Complex_I * std::imag(sigma),
                      reinterpret_cast<_Complex double*>(workev),
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex double*>(resid), ncv,
                      reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex double*>(workd),
                      reinterpret_cast<_Complex double*>(workl), lworkl,
                      reinterpret_cast<_Complex double*>(rwork), &info);
}
}  // namespace arpack
#endif
