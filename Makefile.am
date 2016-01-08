SUBDIRS = UTIL SRC . TESTS EXAMPLES

if MPI
SUBDIRS += PARPACK
AM_DISTCHECK_CONFIGURE_FLAGS = --enable-mpi
endif

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = README PARPACK_CHANGES CHANGES DOCUMENTS VISUAL_STUDIO \
detect_arpack_bug.m4 CMakeLists.txt

# Pkgconfig directory
pkgconfigdir = $(libdir)/pkgconfig

# Files to install in Pkgconfig directory
pkgconfig_DATA = arpack.pc
