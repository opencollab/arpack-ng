LDADD=$(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

check_PROGRAMS = dnsimp bug_1315_single bug_1315_double bug_1323 bug_58_double bug_79_double_complex
if ICB
check_PROGRAMS += icb_arpack_c
check_PROGRAMS += icb_arpack_cpp
endif

TESTS = $(check_PROGRAMS)

EXTRA_DIST = testA.mtx

dnsimp_SOURCES = dnsimp.f mmio.f debug.h

bug_1315_single_SOURCES = bug_1315_single.c
bug_1315_double_SOURCES = bug_1315_double.c
bug_1315_single_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
bug_1315_double_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la

bug_1323_SOURCES = bug_1323.f

bug_58_double_SOURCES = bug_58_double.f

bug_79_double_complex_SOURCES = bug_79_double_complex.f

if ICB
icb_arpack_c_SOURCES = icb_arpack_c.c
icb_arpack_c_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
icb_arpack_c_CFLAGS = -I$(top_builddir)

icb_arpack_cpp_SOURCES = icb_arpack_cpp.cpp
icb_arpack_cpp_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
icb_arpack_cpp_CPPFLAGS = -I$(top_builddir)
endif
