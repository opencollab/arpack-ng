ACLOCAL_AMFLAGS = -I m4

pkgincludedir = $(includedir)/arpack-ng@ITF64SUFFIX@
pkginclude_HEADERS = debug.h stat.h debugF90.h statF90.h arpackdef.h

AM_DISTCHECK_CONFIGURE_FLAGS =
if MPI
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-mpi
endif
if ICB
AM_DISTCHECK_CONFIGURE_FLAGS += --enable-icb
pkginclude_HEADERS += arpackicb.h
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

pkgconfig_DATA = arpack@LIBSUFFIX@@ITF64SUFFIX@.pc parpack@LIBSUFFIX@@ITF64SUFFIX@.pc arpackSolver@LIBSUFFIX@@ITF64SUFFIX@.pc

# Due to the LIBSUFFIX/ITF64SUFFIX, configure doesn't automatically clean this file:
DISTCLEANFILES = arpack@LIBSUFFIX@@ITF64SUFFIX@.pc parpack@LIBSUFFIX@@ITF64SUFFIX@.pc arpackSolver@LIBSUFFIX@@ITF64SUFFIX@.pc
