ACLOCAL_AMFLAGS = -I m4

if MPI
AM_DISTCHECK_CONFIGURE_FLAGS = --enable-mpi
endif

SUBDIRS = UTIL SRC . TESTS EXAMPLES PARPACK

EXTRA_DIST = README PARPACK_CHANGES CHANGES DOCUMENTS VISUAL_STUDIO \
detect_arpack_bug.m4 CMakeLists.txt arpack-ng-config.cmake.in arpack-ng-config-version.cmake.in

pkgconfig_DATA = arpack@LIBSUFFIX@.pc

DISTCLEANFILES = $(pkgconfig_DATA)
