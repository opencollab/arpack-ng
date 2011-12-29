SUBDIRS = UTIL SRC

lib_LTLIBRARIES = libparpack.la
libparpack_la_SOURCES = 
libparpack_la_LDFLAGS = -version-info 2:0
libparpack_la_LIBADD = \
	$(top_builddir)/SRC/libarpacksrc.la \
	$(top_builddir)/UTIL/libarpackutil.la \
	$(top_builddir)/PARPACK/SRC/MPI/libparpacksrcmpi.la \
	$(top_builddir)/PARPACK/UTIL/MPI/libparpackutilmpi.la \
	$(BLAS_LIBS) $(LAPACK_LIBS) $(FLIBS) $(MPILIBS)
