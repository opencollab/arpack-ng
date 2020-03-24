AUTOMAKE_OPTIONS = subdir-objects # Needed as debug/stat* are not in current directory.

F77 = $(MPIF77)
FFLAGS_SAV = @FFLAGS@
FFLAGS =

PSRC = pcontext.f
SSRC = psnaitr.f psnapps.f psnaup2.f psnaupd.f psneigh.f psneupd.f psngets.f \
       pssaitr.f pssapps.f pssaup2.f pssaupd.f psseigt.f psseupd.f pssgets.f \
       psgetv0.f pslarnv.f psnorm2.f

DSRC = pdnaitr.f pdnapps.f pdnaup2.f pdnaupd.f pdneigh.f pdneupd.f pdngets.f \
       pdsaitr.f pdsapps.f pdsaup2.f pdsaupd.f pdseigt.f pdseupd.f pdsgets.f \
       pdgetv0.f pdlarnv.f pdnorm2.f

CSRC = pcnaitr.f pcnapps.f pcnaup2.f pcnaupd.f pcneigh.f pcneupd.f pcngets.f \
       pcgetv0.f pclarnv.f pscnorm2.f

ZSRC = pznaitr.f pznapps.f pznaup2.f pznaupd.f pzneigh.f pzneupd.f pzngets.f \
       pzgetv0.f pzlarnv.f pdznorm2.f

if ICB
SSRC += icbpss.F90 icbpsn.F90
DSRC += icbpds.F90 icbpdn.F90
CSRC += icbpcn.F90
ZSRC += icbpzn.F90
endif

EXTRA_DIST = debug.h stat.h pcontext.h

noinst_LTLIBRARIES = libparpack@LIBSUFFIX@_noopt.la
libparpack@LIBSUFFIX@_noopt_la_SOURCES = pslamch10.f pdlamch10.f
libparpack@LIBSUFFIX@_noopt_la_FFLAGS = -O0

lib_LTLIBRARIES = libparpack@LIBSUFFIX@.la
libparpack@LIBSUFFIX@_la_SOURCES = $(PSRC) $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
libparpack@LIBSUFFIX@_la_SOURCES += $(top_builddir)/dbgini.f
libparpack@LIBSUFFIX@_la_SOURCES += $(top_builddir)/staini.f
libparpack@LIBSUFFIX@_la_FFLAGS = $(FFLAGS_SAV)
libparpack@LIBSUFFIX@_la_LIBADD = libparpack@LIBSUFFIX@_noopt.la \
	$(top_builddir)/PARPACK/UTIL/MPI/libparpackutil.la \
	$(top_builddir)/SRC/libarpack@LIBSUFFIX@.la \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(MPI_Fortran_LIBS)
libparpack@LIBSUFFIX@_la_LDFLAGS = -no-undefined -version-info 3:0:1
if ICB
libparpack@LIBSUFFIX@_la_LIBADD += $(top_builddir)/ICB/libdbgicb.la $(top_builddir)/ICB/libstaicb.la
libparpack@LIBSUFFIX@_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)
endif

pkgconfig_DATA = parpack@LIBSUFFIX@.pc

# Due to the LIBSUFFIX, configure doesn't automatically clean this file:
DISTCLEANFILES = parpack@LIBSUFFIX@.pc
