SUBDIRS = UTIL SRC @MPIDIR@
lib_LTLIBRARIES = libarpack.la
ACLOCAL_AMFLAGS = -I m4/
libarpack_la_SOURCES =
libarpack_la_LDFLAGS = -version-info 2:0
libarpack_la_LIBADD = \
	$(top_builddir)/SRC/libarpacksrc.la \
	$(top_builddir)/UTIL/libarpackutil.la \
	$(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS)
EXTRA_DIST = README


