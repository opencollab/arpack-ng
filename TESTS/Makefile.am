LDADD=$(top_builddir)/libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

AM_DEFAULT_SOURCE_EXT = .f

check_PROGRAMS = dnsimp bug_1315_single bug_1315_double
TESTS = $(check_PROGRAMS)

EXTRA_DIST = testA.mtx

dnsimp_SOURCES = dnsimp.f mmio.f debug.h

bug_1315_single_SOURCES = bug_1315_single.c
bug_1315_double_SOURCES = bug_1315_double.c
