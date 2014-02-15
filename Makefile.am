SUBDIRS = UTIL SRC . TESTS EXAMPLES

if MPI
SUBDIRS += PARPACK
endif

lib_LTLIBRARIES = libarpack.la

ACLOCAL_AMFLAGS = -I m4/

libarpack_la_SOURCES =

# Force libarpack to be linked with the Fortran compiler.
# The file dummy.f does not need to exist in the source tree.
nodist_EXTRA_libarpack_la_SOURCES = dummy.f

libarpack_la_LDFLAGS = -no-undefined -version-info 2:0

libarpack_la_LIBADD = \
	$(top_builddir)/SRC/libarpacksrc.la \
	$(top_builddir)/UTIL/libarpackutil.la \
	$(BLAS_LIBS) $(LAPACK_LIBS)

EXTRA_DIST = README PARPACK_CHANGES CHANGES DOCUMENTS EXAMPLES VISUAL_STUDIO \
detect_arpack_bug.m4

# Pkgconfig directory
pkgconfigdir = $(libdir)/pkgconfig

# Files to install in Pkgconfig directory
pkgconfig_DATA = arpack.pc

