LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

SVD = ssvd dsvd

check_PROGRAMS = $(SVD)
TESTS = $(check_PROGRAMS)

EXTRA_DIST = README debug.h

ssvd_SOURCES = ssvd.f
dsvd_SOURCES = dsvd.f
