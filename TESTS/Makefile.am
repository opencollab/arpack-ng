TESTS = dnsimp bug_1315_single bug_1315_double

check_PROGRAMS = $(TESTS)

AM_DEFAULT_SOURCE_EXT = .f

LDADD=../libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

dnsimp_SOURCES = dnsimp.f mmio.f debug.h

bug_1315_single_SOURCES = bug_1315_single.c

bug_1315_double_SOURCES = bug_1315_double.c

EXTRA_DIST = testA.mtx

