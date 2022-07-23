AUTOMAKE_OPTIONS = subdir-objects # Needed as debug/stat* are not in current directory.

F77 = $(MPIF77)
FFLAGS_SAV = @FFLAGS@
FFLAGS =

PSRC = pcontext.F90
SSRC = psnaitr.F90 psnapps.F90 psnaup2.F90 psnaupd.F90 psneigh.F90 psneupd.F90 psngets.F90 \
       pssaitr.f pssapps.f pssaup2.f pssaupd.f psseigt.f psseupd.f pssgets.f \
       psgetv0.F90 pslarnv.F90 psnorm2.f

DSRC = pdnaitr.F90 pdnapps.F90 pdnaup2.F90 pdnaupd.F90 pdneigh.F90 pdneupd.F90 pdngets.F90 \
       pdsaitr.F90 pdsapps.F90 pdsaup2.F90 pdsaupd.F90 pdseigt.F90 pdseupd.F90 pdsgets.F90 \
       pdgetv0.F90 pdlarnv.F90 pdnorm2.F90

CSRC = pcnaitr.F90 pcnapps.F90 pcnaup2.F90 pcnaupd.F90 pcneigh.F90 pcneupd.F90 pcngets.F90 \
       pcgetv0.F90 pclarnv.F90 pscnorm2.F90

ZSRC = pznaitr.f pznapps.f pznaup2.f pznaupd.f pzneigh.f pzneupd.f pzngets.f \
       pzgetv0.f pzlarnv.f pdznorm2.F90

if ICB
SSRC += icbpss.F90 icbpsn.F90
DSRC += icbpds.F90 icbpdn.F90
CSRC += icbpcn.F90
ZSRC += icbpzn.F90
endif

noinst_LTLIBRARIES = libparpack@LIBSUFFIX@@ITF64SUFFIX@_noopt.la
libparpack@LIBSUFFIX@@ITF64SUFFIX@_noopt_la_SOURCES = pslamch10.F90 pdlamch10.F90
libparpack@LIBSUFFIX@@ITF64SUFFIX@_noopt_la_FFLAGS = -O0

lib_LTLIBRARIES = libparpack@LIBSUFFIX@@ITF64SUFFIX@.la
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_SOURCES = $(PSRC) $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_SOURCES += $(top_builddir)/dbgini.f
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_SOURCES += $(top_builddir)/staini.f
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_FFLAGS = $(FFLAGS_SAV) -I$(top_builddir)
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_FCFLAGS = $(FFLAGS_SAV) -I$(top_builddir)
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_LIBADD = libparpack@LIBSUFFIX@@ITF64SUFFIX@_noopt.la \
	$(top_builddir)/PARPACK/UTIL/MPI/libparpackutil.la \
	$(top_builddir)/SRC/libarpack@LIBSUFFIX@@ITF64SUFFIX@.la \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(MPI_Fortran_LIBS)
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_LDFLAGS = -no-undefined -version-info 3:0:1
if ICB
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_LIBADD += $(top_builddir)/ICB/libdbgicb.la $(top_builddir)/ICB/libstaicb.la
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)
endif
