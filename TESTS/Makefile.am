LDADD=$(top_builddir)/libarpack.la $(LAPACK_LIBS) $(BLAS_LIBS)

check_PROGRAMS = dnsimp bug_1315_single bug_1315_double bug_1323
TESTS = test-dnsimp.sh bug_1315_single bug_1315_double bug_1323

EXTRA_DIST = testA.mtx test-dnsimp.sh

dnsimp_SOURCES = dnsimp.f mmio.f debug.h

bug_1315_single_SOURCES = bug_1315_single.c
bug_1315_double_SOURCES = bug_1315_double.c
bug_1315_single_LDADD = $(top_builddir)/libarpack.la
bug_1315_double_LDADD = $(top_builddir)/libarpack.la

bug_1323_SOURCES = bug_1323.f
