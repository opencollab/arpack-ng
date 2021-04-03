AC_PREREQ(2.67)
AC_INIT([ARPACK-NG],[3.9.0],[https://github.com/opencollab/arpack-ng/issues/],[arpack-ng],[https://github.com/opencollab/arpack-ng/])
AC_CONFIG_SRCDIR([SRC/version.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE

dnl Checks for standard programs.
dnl Autoconf was designed with the fundamental fact that compiler (= language support) is NOT an option.
AC_PROG_FC
AC_PROG_CC
AC_PROG_CXX
PKG_PROG_PKG_CONFIG

dnl Check for BLAS libraries
AX_BLAS([], [AC_MSG_ERROR([cannot find BLAS libraries])])

dnl Check for LAPACK libraries
AX_LAPACK([], [AC_MSG_ERROR([cannot find LAPACK libraries])])

AC_ARG_VAR(INTERFACE64, [set to 1 to use the 64-bit integer interface (ILP64) for ARPACK, BLAS and LAPACK])
if test x"$INTERFACE64" == x"1"; then
  AC_LANG_PUSH([Fortran 77])
  AX_CHECK_COMPILE_FLAG(-fdefault-integer-8, FFLAGS="$FFLAGS -fdefault-integer-8",
                        AX_CHECK_COMPILE_FLAG(-i8, FFLAGS="$FFLAGS -i8",
                                              AC_MSG_WARN([configure does not know how to make your Fortran compiler use 64-bit integers: set it manually via FFLAGS.])))
  AC_LANG_POP([Fortran 77])

  AC_LANG_PUSH([Fortran])
  AX_CHECK_COMPILE_FLAG(-fdefault-integer-8, FCFLAGS="$FCFLAGS -fdefault-integer-8",
                        AX_CHECK_COMPILE_FLAG(-i8, FCFLAGS="$FCFLAGS -i8",
                                              AC_MSG_WARN([configure does not know how to make your Fortran compiler use 64-bit integers: set it manually via FCFLAGS.])))
  AC_LANG_POP([Fortran])
else
  INTERFACE64=0
fi

AC_ARG_VAR(LIBSUFFIX, [suffix to add to ARPACK libraries names])
AC_ARG_VAR(SYMBOLSUFFIX, [suffix to add to ARPACK, BLAS and LAPACK function names])

ifdef([LT_INIT], [], [
  errprint([error: you must have libtool 2.4.2 or a more recent version
])
  m4exit([1])])

LT_PREREQ([2.4.2])
LT_INIT([win32-dll disable-static])

dnl Reminder: checking for ISO_C_BINDING must be done before MPI as MPI changes compilers.

dnl See if compiling with ISO_C_BINDING support
AC_ARG_ENABLE([icb],
              [AS_HELP_STRING([--enable-icb],
                              [build support for *[ae]upd_c with ISO_C_BINDING])],
              [],
              [AS_VAR_SET([enable_icb], [no])])

if test x"$enable_icb" != x"no"; then
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([for iso_c_binding module compilation])
  AC_LINK_IFELSE(
                 [
                  AC_LANG_PROGRAM(
                                  [],
                                  [
                                   USE iso_c_binding
                                   IMPLICIT NONE
                                   INTEGER(C_INT) :: a
                                   a = 1
                                  ]
                                 )
                 ],
                 [AC_MSG_RESULT([yes])],
                 [
                  AC_MSG_RESULT([no])
                  AC_MSG_ERROR([Fortran compiler does not support iso_c_binding])
                 ]
                )
  AC_LANG_POP([Fortran])
fi

dnl See if compiling matrix market example based on ICB
AC_ARG_ENABLE([icb-exmm],
              [AS_HELP_STRING([--enable-icb-exmm],
                              [matrix market example based on ISO_C_BINDING])],
              [AS_VAR_SET([enable_icb], [$enableval])],
              [AS_VAR_SET([enable_icb_exmm], [no])])

if test x"$enable_icb_exmm" != x"no"; then
	PKG_CHECK_MODULES([EIGEN3], [eigen3 >= 3.3])
	AC_LANG_PUSH([C++])
	CPPFLAGS_SAVE=$CPPFLAGS
	CPPFLAGS=$EIGEN3_CFLAGS
	AC_CHECK_HEADER([Eigen/Sparse],
		[],
		[AC_MSG_ERROR([Eigen/Sparse not found])])
	AC_CHECK_HEADER([Eigen/IterativeLinearSolvers],
		[],
		[AC_MSG_ERROR([Eigen/IterativeLinearSolvers not found])])
	AC_CHECK_HEADER([Eigen/SparseLU],
		[],
		[AC_MSG_ERROR([Eigen/SparseLU not found])])
	AC_CHECK_HEADER([Eigen/SparseQR],
		[],
		[AC_MSG_ERROR([Eigen/SparseQR not found])])
	AC_CHECK_HEADER([Eigen/SparseCholesky],
		[],
		[AC_MSG_ERROR([Eigen/SparseCholesky not found])])
	AC_CHECK_HEADER([Eigen/Dense],
		[],
		[AC_MSG_ERROR([Eigen/Dense not found])])
	AC_CHECK_HEADER([Eigen/LU],
		[],
		[AC_MSG_ERROR([Eigen/LU not found])])
	AC_CHECK_HEADER([Eigen/QR],
		[],
		[AC_MSG_ERROR([Eigen/QR not found])])
	AC_CHECK_HEADER([Eigen/Cholesky],
		[],
		[AC_MSG_ERROR([Eigen/Cholesky not found])])
	CPPFLAGS=$CPPFLAGS_SAVE
	AC_LANG_POP([C++])
fi

dnl See if compiling parpack
AC_ARG_ENABLE([mpi],
	[AS_HELP_STRING(
		[--enable-mpi],
		[build parallel version of arpack with MPI])],
	[],
	[AS_VAR_SET([enable_mpi], [no])])
AS_IF([test x"$enable_mpi" != x"no"], [
	FCFLAGS_SAVE=$FCFLAGS
	FCFLAGS=""
	FFLAGS_SAVE=$FFLAGS
	FFLAGS=""

	AC_LANG_PUSH([Fortran 77])
	AX_MPI([], AC_MSG_ERROR([could not compile a Fortran MPI test program]))
	AC_SUBST([MPI_Fortran_LIBS], [$MPILIBS])
	F77=$MPIF77
	AC_LANG_POP([Fortran 77])

	FCFLAGS=$FCFLAGS_SAVE
	FFLAGS=$FFLAGS_SAVE
])

if test x"$enable_icb" != x"no"; then
  if test x"$enable_mpi" != x"no"; then
    FCFLAGS_SAVE=$FCFLAGS
    FCFLAGS=""
    FFLAGS_SAVE=$FFLAGS
    FFLAGS=""

    AC_LANG_PUSH([Fortran])
    AX_MPI([], AC_MSG_ERROR([could not compile a Fortran MPI test program]))
    FC=$MPIFC
    AC_FC_LIBRARY_LDFLAGS dnl set FCLIBS
    AC_SUBST([MPI_Fortran_LIBS], ["$MPILIBS $FCLIBS"])
    AC_LANG_POP([Fortran])

    AC_LANG_PUSH([C])
    AX_MPI([], AC_MSG_ERROR([could not compile a C MPI test program]))
    AC_SUBST([MPI_C_LIBS], ["$MPILIBS $FCLIBS"])
    CC=$MPICC
    AC_LANG_POP([C])

    AC_LANG_PUSH([C++])
    AX_MPI([], AC_MSG_ERROR([could not compile a C++ MPI test program]))
    AC_SUBST([MPI_CXX_LIBS], ["$MPILIBS $FCLIBS"])
    CXX=$MPICXX
    AC_LANG_POP([C++])

    AC_LANG_PUSH([C])
    AC_MSG_CHECKING([for MPI_Comm_c2f function])
    AC_LINK_IFELSE(
                   [
                    AC_LANG_PROGRAM(
                                    [[#include <mpi.h>]],
                                    [[MPI_Comm c_comm; MPI_Comm_c2f(c_comm)]]
                                   )
                   ],
                   [AC_MSG_RESULT([yes])],
                   [
                    AC_MSG_RESULT([no])
                    AC_MSG_ERROR([symbol MPI_Comm_c2f does not exist])
                   ]
                  )
    AC_LANG_POP([C])

    FCFLAGS=$FCFLAGS_SAVE
    FFLAGS=$FFLAGS_SAVE
  fi
  AX_CXX_COMPILE_STDCXX(11)
else
  # ICB saves you from old-fashion-boring-cumbersome-fortran/C crap... For arpack symbols (only).

  if test x"$SYMBOLSUFFIX" != x""; then
    dnl Need to rely on non-F77 features
    AC_FC_LINE_LENGTH(unlimited)
    FFLAGS="$FFLAGS $ac_cv_fc_line_length"

    AC_LANG_PUSH([Fortran 77])
    AX_CHECK_COMPILE_FLAG(-cpp, FFLAGS="$FFLAGS -cpp",
                          AX_CHECK_COMPILE_FLAG(-fpp, FFLAGS="$FFLAGS -fpp",
                                                AC_MSG_WARN([configure does not know how to enable your Fortran compiler's preprocessor: set it manually via FFLAGS.])))
    AC_LANG_POP([Fortran 77])

    SCALARFUNS="naitr napps naup2 naupd nconv neigh neupd ngets statn saitr sapps saup2 saupd sconv seigt seupd sgets stats getv0 sortc sortr sesrt stqrb"
    COMPLEXFUNS="naitr napps naup2 naupd neigh neupd ngets statn getv0 sortc"

    FUNS1="$SCALARFUNS"
    for f in $FUNS1
    do
      FFLAGS="$FFLAGS -Ds$f=s$f$SYMBOLSUFFIX -Dd$f=d$f$SYMBOLSUFFIX"
    done

    FUNS2="$COMPLEXFUNS"
    for f in $FUNS2
    do
      FFLAGS="$FFLAGS -Dc$f=c$f$SYMBOLSUFFIX -Dz$f=z$f$SYMBOLSUFFIX"
    done
  fi

  dnl Both a lowercase and uppercase variable are needed
  dnl since the macro changes the case of the whole expression
  symbolsuffix=$SYMBOLSUFFIX
  AC_F77_FUNC(snaupd$SYMBOLSUFFIX, snaupdsuff)
  AC_F77_FUNC(sneupd$SYMBOLSUFFIX, sneupdsuff)
  AC_F77_FUNC(dnaupd$SYMBOLSUFFIX, dnaupdsuff)
  AC_F77_FUNC(dneupd$SYMBOLSUFFIX, dneupdsuff)
  dnl Needed for BLAS check and for tests (even when suffix is empty)
  CPPFLAGS="$CPPFLAGS -Dsnaupd=$snaupdsuff -Dsneupd=$sneupdsuff -Ddnaupd=$dnaupdsuff -Ddneupd=$dneupdsuff"
fi

dnl As BLAS/LAPACK does not provide ICB, we may have to deal with symbols the old-fashion-boring-cumbersome-fortran/C way...

if test x"$SYMBOLSUFFIX" != x""; then
  dnl Need to rely on non-F77 features
  AC_FC_LINE_LENGTH(unlimited)
  FFLAGS="$FFLAGS $ac_cv_fc_line_length"

  AC_LANG_PUSH([Fortran 77])
  AX_CHECK_COMPILE_FLAG(-cpp, FFLAGS="$FFLAGS -cpp",
                        AX_CHECK_COMPILE_FLAG(-fpp, FFLAGS="$FFLAGS -fpp",
                                              AC_MSG_WARN([configure does not know how to enable your Fortran compiler's preprocessor: set it manually via FFLAGS.])))
  AC_LANG_POP([Fortran 77])

  BLASFUNS1="axpy copy gemv geqr2 lacpy lae2 lahqr lanhs larnv lartg lascl laset lasrt scal trevc trmm trsen gbmv gbtrf gbtrs gttrf gttrs pttrf pttrs"
  BLASFUNS2="dot ger labad laev2 lamch lanst lanv2 lapy2 larf larfg lasr nrm2 orm2r rot steqr swap"
  BLASFUNS3="dotc geru unm2r"
  BLASFUNS4="COPY LABAD LAMCH LANHS LANV2 LARFG ROT GEMV"
  BLASFUNS5="scnrm2 dznrm2 csscal zdscal"

  FUNS1="$BLASFUNS1 $BLASFUNS2"
  for f in $FUNS1
  do
    FFLAGS="$FFLAGS -Ds$f=s$f$SYMBOLSUFFIX -Dd$f=d$f$SYMBOLSUFFIX"
  done

  FUNS2="$BLASFUNS1 $BLASFUNS3"
  for f in $FUNS2
  do
    FFLAGS="$FFLAGS -Dc$f=c$f$SYMBOLSUFFIX -Dz$f=z$f$SYMBOLSUFFIX"
  done

  for f in $BLASFUNS4
  do
    FFLAGS="$FFLAGS -DS$f=S$f$SYMBOLSUFFIX -DD$f=D$f$SYMBOLSUFFIX"
  done

  for f in $BLASFUNS5
  do
    FFLAGS="$FFLAGS -D$f=$f$SYMBOLSUFFIX -D$f=$f$SYMBOLSUFFIX"
  done

  dnl Both a lowercase and uppercase variable are needed
  dnl since the macro changes the case of the whole expression
  symbolsuffix=$SYMBOLSUFFIX
  AC_F77_FUNC(sgemm$SYMBOLSUFFIX, sgemmsuff)
  AC_F77_FUNC(cheev$SYMBOLSUFFIX, cheevsuff)
  dnl Needed for BLAS check and for tests (even when suffix is empty)
  CPPFLAGS="$CPPFLAGS -Dsgemm=$sgemmsuff -Dcheev=$cheevsuff"
fi

dnl Set conditional for Makefiles.

AM_CONDITIONAL([MPI], [test x"$enable_mpi" != x"no"])
AM_CONDITIONAL([BLACS], [false])
AM_CONDITIONAL([ICB], [test x"$enable_icb" != x"no"])
AM_CONDITIONAL([ICBEXMM], [test x"$enable_icb_exmm" != x"no"])

m4_ifdef([PKG_INSTALLDIR], [PKG_INSTALLDIR], [AC_SUBST([pkgconfigdir], [${libdir}/pkgconfig])])
AC_SUBST([ARPACK_PC_LIBS_PRIVATE], ["$LAPACK_LIBS $BLAS_LIBS"])
AC_SUBST([PARPACK_PC_LIBS_PRIVATE], ["$LAPACK_LIBS $BLAS_LIBS $MPI_Fortran_LIBS"])
AC_CONFIG_FILES([
	SRC/arpack$LIBSUFFIX.pc:SRC/arpack.pc.in
	PARPACK/SRC/MPI/parpack$LIBSUFFIX.pc:PARPACK/SRC/MPI/parpack.pc.in
], [], [LIBSUFFIX="$LIBSUFFIX"])

AC_CONFIG_FILES([
	EXAMPLES/MATRIX_MARKET/arpackSolver.pc
	arpackdef.h
	arpackicb.h
	arpack-ng-config-version.cmake
	arpack-ng-config.cmake
	Makefile
	ICB/Makefile
	UTIL/Makefile
	SRC/Makefile
	TESTS/Makefile
	EXAMPLES/Makefile
	EXAMPLES/BAND/Makefile
	EXAMPLES/COMPLEX/Makefile
	EXAMPLES/NONSYM/Makefile
	EXAMPLES/SIMPLE/Makefile
	EXAMPLES/SVD/Makefile
	EXAMPLES/SYM/Makefile
	EXAMPLES/MATRIX_MARKET/Makefile
	PARPACK/Makefile
	PARPACK/SRC/MPI/Makefile
	PARPACK/SRC/BLACS/Makefile
	PARPACK/UTIL/MPI/Makefile
	PARPACK/UTIL/BLACS/Makefile
	PARPACK/TESTS/MPI/Makefile
	PARPACK/EXAMPLES/MPI/Makefile
	PARPACK/EXAMPLES/BLACS/Makefile
])
AC_CONFIG_LINKS([
	TESTS/testA.mtx:TESTS/testA.mtx
	EXAMPLES/MATRIX_MARKET/An.mtx:EXAMPLES/MATRIX_MARKET/An.mtx
	EXAMPLES/MATRIX_MARKET/As.mtx:EXAMPLES/MATRIX_MARKET/As.mtx
	EXAMPLES/MATRIX_MARKET/Az.mtx:EXAMPLES/MATRIX_MARKET/Az.mtx
	EXAMPLES/MATRIX_MARKET/B.mtx:EXAMPLES/MATRIX_MARKET/B.mtx
	EXAMPLES/MATRIX_MARKET/Bz.mtx:EXAMPLES/MATRIX_MARKET/Bz.mtx
])
AC_OUTPUT

AC_MSG_RESULT([
--------------------------------------------------
Configuration summary for $PACKAGE_STRING
--------------------------------------------------
Installation prefix : $prefix
MPI enabled         : $enable_mpi
ICB enabled         : $enable_icb
INTERFACE64         : $INTERFACE64
F77                 : $F77
FFLAGS              : $FFLAGS
FC                  : $FC
FCFLAGS             : $FCFLAGS
CC                  : $CC
CFLAGS              : $CFLAGS
CXX                 : $CXX
CXXFLAGS            : $CXXFLAGS
CPPFLAGS            : $CPPFLAGS
MPI_Fortran_LIBS    : $MPI_Fortran_LIBS
MPI_C_LIBS          : $MPI_C_LIBS
MPI_CXX_LIBS        : $MPI_CXX_LIBS
BLAS                : $BLAS_LIBS
LAPACK              : $LAPACK_LIBS
EIGEN               : $EIGEN3_CFLAGS
LIBS                : $LIBS
LDADD               : $LDADD
--------------------------------------------------
Configuration OK
--------------------------------------------------
])
