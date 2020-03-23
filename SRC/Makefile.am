AUTOMAKE_OPTIONS = subdir-objects # Needed as debug/stat* are not in current directory.

SSRC = snaitr.f snapps.f snaup2.f snaupd.f snconv.f sneigh.f sneupd.f sngets.f sstatn.f \
       ssaitr.f ssapps.f ssaup2.f ssaupd.f ssconv.f sseigt.f sseupd.f ssgets.f sstats.f \
       sgetv0.f ssortc.f ssortr.f ssesrt.f sstqrb.f

DSRC = dnaitr.f dnapps.f dnaup2.f dnaupd.f dnconv.f dneigh.f dneupd.f dngets.f dstatn.f \
       dsaitr.f dsapps.f dsaup2.f dsaupd.f dsconv.f dseigt.f dseupd.f dsgets.f dstats.f \
       dgetv0.f dsortc.f dsortr.f dsesrt.f dstqrb.f

CSRC = cnaitr.f cnapps.f cnaup2.f cnaupd.f cneigh.f cneupd.f cngets.f cstatn.f \
       cgetv0.f csortc.f

ZSRC = znaitr.f znapps.f znaup2.f znaupd.f zneigh.f zneupd.f zngets.f zstatn.f \
       zgetv0.f zsortc.f

if ICB
SSRC += icbass.F90 icbasn.F90
DSRC += icbads.F90 icbadn.F90
CSRC += icbacn.F90
ZSRC += icbazn.F90
endif

EXTRA_DIST = debug.h stat.h version.h

lib_LTLIBRARIES = libarpack@LIBSUFFIX@.la
libarpack@LIBSUFFIX@_la_SOURCES = $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
libarpack@LIBSUFFIX@_la_SOURCES += $(top_builddir)/dbgini.f
libarpack@LIBSUFFIX@_la_SOURCES += $(top_builddir)/staini.f
libarpack@LIBSUFFIX@_la_LIBADD = $(top_builddir)/UTIL/libarpackutil.la $(LAPACK_LIBS) $(BLAS_LIBS)
libarpack@LIBSUFFIX@_la_LDFLAGS = -no-undefined -version-info 3:0:1
if ICB
libarpack@LIBSUFFIX@_la_LIBADD += $(top_builddir)/ICB/libdbgicb.la $(top_builddir)/ICB/libstaicb.la
libarpack@LIBSUFFIX@_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir)
endif

pkgconfig_DATA = arpack@LIBSUFFIX@.pc

# Due to the LIBSUFFIX, configure doesn't automatically clean this file:
DISTCLEANFILES = arpack@LIBSUFFIX@.pc
