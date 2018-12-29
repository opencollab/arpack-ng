ACLOCAL_AMFLAGS = -I m4

arpack_includedir = $(includedir)/arpack

AM_DISTCHECK_CONFIGURE_FLAGS =
if MPI
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-mpi
endif
arpack_include_HEADERS = debug.h stat.h arpackdef.h

SUBDIRS = .
if ICB
SUBDIRS += ICB
endif
SUBDIRS += UTIL SRC TESTS EXAMPLES
if MPI
SUBDIRS += PARPACK
endif

DIST_SUBDIRS = .
if ICB
DIST_SUBDIRS += ICB
endif
DIST_SUBDIRS += UTIL SRC TESTS EXAMPLES
if MPI
DIST_SUBDIRS += PARPACK
endif

EXTRA_DIST = README.md PARPACK_CHANGES CHANGES DOCUMENTS VISUAL_STUDIO \
detect_arpack_bug.m4 CMakeLists.txt arpack-ng-config.cmake.in arpack-ng-config-version.cmake.in

pkgconfig_DATA = arpack@LIBSUFFIX@.pc

DISTCLEANFILES = $(pkgconfig_DATA) arpackdef.h
