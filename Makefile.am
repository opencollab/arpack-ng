ACLOCAL_AMFLAGS = -I m4

pkgincludedir = $(includedir)/arpack
pkginclude_HEADERS = debug.h stat.h arpackdef.h

AM_DISTCHECK_CONFIGURE_FLAGS =
if MPI
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-mpi
endif
if ICB
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-icb
endif
if ICBEXMM
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-icb-exmm
endif

SUBDIRS = .
if ICB
SUBDIRS += ICB
endif
SUBDIRS += UTIL SRC TESTS EXAMPLES
if MPI
SUBDIRS += PARPACK
endif

EXTRA_DIST = README.md PARPACK_CHANGES CHANGES DOCUMENTS VISUAL_STUDIO \
detect_arpack_bug.m4 CMakeLists.txt

cmakedir = $(libdir)/cmake/arpack-ng
cmake_DATA = arpack-ng-config.cmake arpack-ng-config-version.cmake
