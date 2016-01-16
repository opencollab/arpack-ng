AC_PREREQ(2.59)
AC_INIT([ARPACK-NG],[3.4.0],[https://github.com/opencollab/arpack-ng/issues/],[arpack-ng],[https://github.com/opencollab/arpack-ng/])
AM_INIT_AUTOMAKE([foreign])
AM_MAINTAINER_MODE

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_LINKS([TESTS/testA.mtx:TESTS/testA.mtx])

dnl Checks for standard programs.
AC_PROG_F77
AC_PROG_CC

ifdef([LT_INIT], [], [
  errprint([error: you must have libtool 2.4.2 or a more recent version
])
  m4exit([1])])

LT_PREREQ([2.4.2])
LT_INIT([win32-dll])

AC_ARG_VAR(LIBSUFFIX, [suffix to add to ARPACK library name])
AC_ARG_VAR(SYMBOLSUFFIX, [suffix to add to ARPACK, BLAS and LAPACK function names])

if test x"$SYMBOLSUFFIX" != x""; then
  FFLAGS+=" -cpp -ffixed-line-length-none"

  SCALARFUNS="naitr napps naup2 naupd nconv neigh neupd ngets statn saitr sapps saup2 saupd sconv seigt seupd sgets stats getv0 sortc sortr sesrt stqrb"
  COMPLEXFUNS="naitr napps naup2 naupd neigh neupd ngets statn getv0 sortc"

  BLASFUNS1="axpy copy gemv geqr2 lacpy lae2 lahqr lanhs larnv lartg lascl laset lasrt scal trevc trmm trsen gbmv gbtrf gbtrs gttrf gttrs pttrf pttrs"
  BLASFUNS2="dot ger labad laev2 lamch lanst lanv2 lapy2 larf larfg lasr nrm2 orm2r rot steqr swap"
  BLASFUNS3="dotc geru unm2r"
  BLASFUNS4="COPY LABAD LAMCH LANHS LANV2 LARFG ROT GEMV"

  FUNS1="$SCALARFUNS $BLASFUNS1 $BLASFUNS2"
  for f in $FUNS1
  do
    FFLAGS+=" -Ds$f=s$f$SYMBOLSUFFIX -Dd$f=d$f$SYMBOLSUFFIX"
  done

  FUNS2="$COMPLEXFUNS $BLASFUNS1 $BLASFUNS3"
  for f in $FUNS2
  do
    FFLAGS+=" -Dc$f=c$f$SYMBOLSUFFIX -Dz$f=z$f$SYMBOLSUFFIX"
  done

  for f in $BLASFUNS4
  do
    FFLAGS+=" -DS$f=S$f$SYMBOLSUFFIX -DD$f=D$f$SYMBOLSUFFIX"
  done

  FFLAGS+=" -Dscnrm2=scnrm2$SYMBOLSUFFIX -Ddznrm2=dznrm2$SYMBOLSUFFIX -Dcsscal=csscal$SYMBOLSUFFIX -Dzdscal=zdscal$SYMBOLSUFFIX"

  dnl Needed for BLAS check and for tests
  CFLAGS+=" -Dsgemm_=sgemm${SYMBOLSUFFIX}_ -Dsnaupd_=snaupd${SYMBOLSUFFIX}_ -Dsneupd_=sneupd${SYMBOLSUFFIX}_ -Ddnaupd_=dnaupd${SYMBOLSUFFIX}_ -Ddneupd_=dneupd${SYMBOLSUFFIX}_ -Dcheev_=cheev${SYMBOLSUFFIX}_"
fi

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
