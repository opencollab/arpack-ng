AC_PREREQ(2.67)
AC_INIT([ARPACK-NG],[3.5.0],[https://github.com/opencollab/arpack-ng/issues/],[arpack-ng],[https://github.com/opencollab/arpack-ng/])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_LINKS([TESTS/testA.mtx:TESTS/testA.mtx])

dnl Checks for standard programs.
AC_PROG_F77
AC_PROG_CC

dnl Check for BLAS libraries
sinclude(ax_blas.m4)
AX_BLAS
if test "$ax_blas_ok" = "no"; then
   AC_MSG_ERROR([cannot find BLAS libraries])
fi

dnl Check for LAPACK libraries
sinclude(ax_lapack.m4)
AX_LAPACK
if test "$ax_lapack_ok" = "no"; then
   AC_MSG_ERROR([cannot find LAPACK libraries])
fi

AC_ARG_VAR(INTERFACE64, [set to 1 to use the 64-bit integer interface (ILP64) for ARPACK, BLAS and LAPACK])
if test x"$INTERFACE64" == x"1"; then
  AX_CHECK_COMPILE_FLAG(-fdefault-integer-8, FFLAGS="$FFLAGS -fdefault-integer-8",
                        AX_CHECK_COMPILE_FLAG(-i8, FFLAGS="$FFLAGS -i8",
                                              AC_MSG_WARN([configure does not know how to make your Fortran compiler use 64-bit integers: set it manually via FFLAGS.])))
  CFLAGS="$CFLAGS -DBLASINT=int64_t"
fi

AC_ARG_VAR(LIBSUFFIX, [suffix to add to ARPACK libraries names])
AC_ARG_VAR(SYMBOLSUFFIX, [suffix to add to ARPACK, BLAS and LAPACK function names])

if test x"$SYMBOLSUFFIX" != x""; then
  dnl Need to rely on non-F77 features
  AC_FC_LINE_LENGTH(unlimited)
  FFLAGS="$FFLAGS $ac_cv_fc_line_length"

  AX_CHECK_COMPILE_FLAG(-cpp, FFLAGS="$FFLAGS -cpp",
                        AX_CHECK_COMPILE_FLAG(-fpp, FFLAGS="$FFLAGS -fpp",
                                              AC_MSG_WARN([configure does not know how to enable your Fortran compiler's preprocessor: set it manually via FFLAGS.])))

  SCALARFUNS="naitr napps naup2 naupd nconv neigh neupd ngets statn saitr sapps saup2 saupd sconv seigt seupd sgets stats getv0 sortc sortr sesrt stqrb"
  COMPLEXFUNS="naitr napps naup2 naupd neigh neupd ngets statn getv0 sortc"

  BLASFUNS1="axpy copy gemv geqr2 lacpy lae2 lahqr lanhs larnv lartg lascl laset lasrt scal trevc trmm trsen gbmv gbtrf gbtrs gttrf gttrs pttrf pttrs"
  BLASFUNS2="dot ger labad laev2 lamch lanst lanv2 lapy2 larf larfg lasr nrm2 orm2r rot steqr swap"
  BLASFUNS3="dotc geru unm2r"
  BLASFUNS4="COPY LABAD LAMCH LANHS LANV2 LARFG ROT GEMV"
  BLASFUNS5="scnrm2 dznrm2 csscal zdscal"

  FUNS1="$SCALARFUNS $BLASFUNS1 $BLASFUNS2"
  for f in $FUNS1
  do
    FFLAGS="$FFLAGS -Ds$f=s$f$SYMBOLSUFFIX -Dd$f=d$f$SYMBOLSUFFIX"
  done

  FUNS2="$COMPLEXFUNS $BLASFUNS1 $BLASFUNS3"
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
fi

dnl Both a lowercase and uppercase variable are needed
dnl since the macro changes the case of the whole expression
symbolsuffix=$SYMBOLSUFFIX
AC_F77_FUNC(sgemm$SYMBOLSUFFIX, sgemmsuff)
AC_F77_FUNC(snaupd$SYMBOLSUFFIX, snaupdsuff)
AC_F77_FUNC(sneupd$SYMBOLSUFFIX, sneupdsuff)
AC_F77_FUNC(dnaupd$SYMBOLSUFFIX, dnaupdsuff)
AC_F77_FUNC(dneupd$SYMBOLSUFFIX, dneupdsuff)
AC_F77_FUNC(cheev$SYMBOLSUFFIX, cheevsuff)
dnl Needed for BLAS check and for tests (even when suffix is empty)
CFLAGS="$CFLAGS -Dsgemm=$sgemmsuff -Dsnaupd=$snaupdsuff -Dsneupd=$sneupdsuff -Ddnaupd=$dnaupdsuff -Ddneupd=$dneupdsuff -Dcheev=$cheevsuff"

ifdef([LT_INIT], [], [
  errprint([error: you must have libtool 2.4.2 or a more recent version
])
  m4exit([1])])

LT_PREREQ([2.4.2])
LT_INIT([win32-dll])

dnl See if compiling parpack
AC_ARG_ENABLE(mpi,
	      AC_HELP_STRING([--enable-mpi],
			     [build parallel version of arpack with MPI]),
	      [enable_mpi=$enableval],
	      [enable_mpi=no])
if test x"$enable_mpi" != x"no"; then
  AC_LANG_PUSH([Fortran 77])
  AX_MPI([], AC_MSG_ERROR([could not compile a MPI test program]))
  AC_LANG_POP([Fortran 77])
fi

AM_CONDITIONAL([MPI],[test x"$enable_mpi" = x"yes"])
AM_CONDITIONAL([BLACS],[false])

AC_CONFIG_FILES([arpack$LIBSUFFIX.pc:arpack.pc.in], [], [LIBSUFFIX="$LIBSUFFIX"])

AC_CONFIG_FILES([
	Makefile
	SRC/Makefile
	UTIL/Makefile
	TESTS/Makefile
	EXAMPLES/Makefile
	EXAMPLES/BAND/Makefile
	EXAMPLES/COMPLEX/Makefile
	EXAMPLES/NONSYM/Makefile
	EXAMPLES/SIMPLE/Makefile
	EXAMPLES/SVD/Makefile
	EXAMPLES/SYM/Makefile
	PARPACK/Makefile
	PARPACK/SRC/MPI/Makefile
	PARPACK/SRC/BLACS/Makefile
	PARPACK/UTIL/MPI/Makefile
	PARPACK/UTIL/BLACS/Makefile
	PARPACK/EXAMPLES/MPI/Makefile
	PARPACK/EXAMPLES/BLACS/Makefile
	])
AC_OUTPUT
