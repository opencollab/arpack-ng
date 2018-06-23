dnl
dnl Check whether ARPACK works (does not crash)
dnl
dnl Using a pure Fortran program doesn't seem to crash when linked
dnl with the buggy ARPACK library but the C++ program does.  Maybe
dnl it is the memory allocation that exposes the bug and using statically
dnl allocated arrays in Fortran does not?
dnl
dnl This code is not used by arpack-ng itself.
dnl This is a macro for applications using arpack to detect that the version
dnl of arpack behave correctly (ie not arpack-ng)
dnl This is the work of Rik <rik@octave.org>
dnl 
dnl This code is released under the same license as arpack
dnl
AC_DEFUN([CHECK_ARPACK_OK], [
  AC_LANG_PUSH(C++)
  AC_CACHE_CHECK([whether the arpack library works],
    [cv_lib_arpack_ok], [
      AC_RUN_IFELSE([AC_LANG_PROGRAM([[
// External functions from ARPACK library
extern "C" int
F77_FUNC (dnaupd, DNAUPD) (int&, const char *, const int&, const char *,
                           int&, const double&, double*, const int&,
                           double*, const int&, int*, int*, double*,
                           double*, const int&, int&, long int, long int);

extern "C" int
F77_FUNC (dneupd, DNEUPD) (const int&, const char *, int*, double*,
                           double*, double*, const int&,
                           const double&, const double&, double*,
                           const char*, const int&, const char *,
                           int&, const double&, double*, const int&,
                           double*, const int&, int*, int*, double*,
                           double*, const int&, int&, long int,
                           long int, long int);

extern "C" int
F77_FUNC (dgemv, DGEMV) (const char *, const int&, const int&,
                         const double&, const double*, const int&,
                         const double*, const int&, const double&,
                         double*, const int&, long int);

#include <cfloat>

void
doit (void)
{
  // Based on the octave function EigsRealNonSymmetricMatrix from
  // liboctave/eigs-base.cc.

  // Problem matrix.  See bug #31479
  int n = 4;
  double *m = new double [n * n];
  m[0] = 1, m[4] = 0, m[8]  = 0, m[12] = -1;
  m[1] = 0, m[5] = 1, m[9]  = 0, m[13] = 0;
  m[2] = 0, m[6] = 0, m[10] = 1, m[14] = 0;
  m[3] = 0, m[7] = 0, m[11] = 2, m[15] = 1;

  double *resid = new double [4];

  resid[0] = 0.960966;
  resid[1] = 0.741195;
  resid[2] = 0.150143;
  resid[3] = 0.868067;

  int *ip = new int [11];

  ip[0] = 1;   // ishift
  ip[1] = 0;   // ip[1] not referenced
  ip[2] = 300; // mxiter, maximum number of iterations
  ip[3] = 1;   // NB blocksize in recurrence
  ip[4] = 0;   // nconv, number of Ritz values that satisfy convergence
  ip[5] = 0;   // ip[5] not referenced
  ip[6] = 1;   // mode
  ip[7] = 0;   // ip[7] to ip[10] are return values
  ip[8] = 0;
  ip[9] = 0;
  ip[10] = 0;
 
  int *ipntr = new int [14];

  int k = 1;
  int p = 3;
  int lwork = 3 * p * (p + 2);

  double *v = new double [n * (p + 1)];
  double *workl = new double [lwork + 1];
  double *workd = new double [3 * n + 1];

  int ido = 0;
  int info = 0;

  double tol = DBL_EPSILON;

  do 
    {
      F77_FUNC (dnaupd, DNAUPD) (ido, "I", n, "LM", k, tol, resid, p,
                                 v, n, ip, ipntr, workd, workl, lwork,
                                 info, 1L, 2L);

      if (ido == -1 || ido == 1 || ido == 2)
        {
          double *x = workd + ipntr[0] - 1;
          double *y = workd + ipntr[1] - 1;

          F77_FUNC (dgemv, DGEMV) ("N", n, n, 1.0, m, n, x, 1, 0.0,
                                   y, 1, 1L);
        }
      else
        {
          if (info < 0)
            {
              return;  // Error
            }

          break;
        }
    } 
  while (1);

  int *sel = new int [p];

  // In Octave, the dimensions of dr and di are k+1, but k+2 avoids segfault
  double *dr = new double [k + 1];
  double *di = new double [k + 1];
  double *workev = new double [3 * p];

  for (int i = 0; i < k + 1; i++)
    dr[i] = di[i] = 0.;

  int rvec = 1;

  double sigmar = 0.0;
  double sigmai = 0.0;

  // In Octave, this is n*(k+1), but k+2 avoids segfault
  double *z = new double [n * (k + 1)];

  F77_FUNC (dneupd, DNEUPD) (rvec, "A", sel, dr, di, z, n, sigmar,
                             sigmai, workev, "I", n, "LM", k, tol,
                             resid, p, v, n, ip, ipntr, workd,
                             workl, lwork, info, 1L, 1L, 2L);
}
]], [[
  for (int i = 0; i < 10; i++)
    doit ();
]])],
  [cv_lib_arpack_ok=yes],
  [cv_lib_arpack_ok=no],
  [cv_lib_arpack_ok=yes])])
  AC_LANG_POP(C++)
  if test "$cv_lib_arpack_ok" = "yes"; then
    $1
  else
    $2
  fi
])
