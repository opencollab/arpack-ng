check_PROGRAMS = dnsimp  bug_1323

dnsimp_SOURCES = dnsimp.f mmio.f debug.h

dnsimp_LDADD=../libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

EXTRA_DIST = testA.mtx

bug_1323_SOURCES =  bug_1323.f

bug_1323_LDADD=../libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

TESTS = dnsimp bug_1323
