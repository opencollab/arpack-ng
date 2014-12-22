LDADD = $(top_builddir)/libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

SVD = ssvd dsvd

check_PROGRAMS = $(SVD)
TESTS = $(check_PROGRAMS)

ssvd_SOURCES = ssvd.f
dsvd_SOURCES = dsvd.f
