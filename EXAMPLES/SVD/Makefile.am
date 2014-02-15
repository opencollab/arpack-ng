check_PROGRAMS = dsvd ssvd

dsvd_SOURCES = dsvd.f

dsvd_LDADD=../../libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)


ssvd_SOURCES = ssvd.f

ssvd_LDADD=../../libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

TESTS = dsvd ssvd
