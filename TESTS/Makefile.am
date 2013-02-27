pkglibdirexec_PROGRAMS = dnsimp 

dnsimp_SOURCES = dnsimp.f mmio.f debug.h

dnsimp_LDADD=../libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

EXTRA_DIST = testA.mtx