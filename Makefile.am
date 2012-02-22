SUBDIRS = UTIL SRC PARPACK . TESTS 


lib_LTLIBRARIES = libarpack.la

ACLOCAL_AMFLAGS = -I m4/

libarpack_la_SOURCES =

libarpack_la_LDFLAGS = -version-info 2:0

libarpack_la_LIBADD = \
	$(top_builddir)/SRC/libarpacksrc.la \
	$(top_builddir)/UTIL/libarpackutil.la \
	$(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS)

EXTRA_DIST = README CHANGES

# Pkgconfig directory
pkgconfigdir = $(libdir)/pkgconfig

# Files to install in Pkgconfig directory
pkgconfig_DATA = arpack.pc
