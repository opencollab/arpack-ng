AUTOMAKE_OPTIONS = subdir-objects # Needed as debug/stat* are not in current directory.

F77 = $(MPIF77)

PSRC = pcontext.f
SSRC = psnaitr.f psnapps.f psnaup2.f psnaupd.f psneigh.f psneupd.f psngets.f \
       pssaitr.f pssapps.f pssaup2.f pssaupd.f psseigt.f psseupd.f pssgets.f \
       psgetv0.f pslarnv.f psnorm2.f pslamch10.f

DSRC = pdnaitr.f pdnapps.f pdnaup2.f pdnaupd.f pdneigh.f pdneupd.f pdngets.f \
       pdsaitr.f pdsapps.f pdsaup2.f pdsaupd.f pdseigt.f pdseupd.f pdsgets.f \
       pdgetv0.f pdlarnv.f pdnorm2.f pdlamch10.f

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

lib_LTLIBRARIES = libparpack@LIBSUFFIX@@ITF64SUFFIX@.la
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_SOURCES = $(PSRC) $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_SOURCES += $(top_builddir)/dbgini.f
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_SOURCES += $(top_builddir)/staini.f
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_LIBADD = \
	$(top_builddir)/PARPACK/UTIL/MPI/libparpackutil.la \
	$(top_builddir)/SRC/libarpack@LIBSUFFIX@@ITF64SUFFIX@.la \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(MPI_Fortran_LIBS)
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_LDFLAGS = -no-undefined -version-info 3:0:1
if ICB
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_LIBADD += $(top_builddir)/ICB/libdbgicb.la $(top_builddir)/ICB/libstaicb.la
libparpack@LIBSUFFIX@@ITF64SUFFIX@_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)
endif
