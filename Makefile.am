ACLOCAL_AMFLAGS = -I m4

AM_DISTCHECK_CONFIGURE_FLAGS =
if MPI
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-mpi
endif
include_HEADERS = debug-arpack.h stat.h
if ICB
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-icb
include_HEADERS += debug_c.h debug_c.hpp
include_HEADERS += stat_c.h stat_c.hpp
include_HEADERS += arpack.h arpack.hpp
if MPI
include_HEADERS += parpack.h parpack.hpp
endif
endif

SUBDIRS = UTIL SRC . TESTS EXAMPLES PARPACK

EXTRA_DIST = README.md PARPACK_CHANGES CHANGES DOCUMENTS VISUAL_STUDIO \
detect_arpack_bug.m4 CMakeLists.txt arpack-ng-config.cmake.in arpack-ng-config-version.cmake.in \
debug-arpack.h stat.h
if ICB
EXTRA_DIST += debug_c.h debug_c.hpp
EXTRA_DIST += stat_c.h stat_c.hpp
EXTRA_DIST += arpack.h arpack.hpp
if MPI
EXTRA_DIST += parpack.h parpack.hpp
endif
endif

pkgconfig_DATA = arpack@LIBSUFFIX@.pc

DISTCLEANFILES = $(pkgconfig_DATA)
