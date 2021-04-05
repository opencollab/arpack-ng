LDADD=$(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

check_PROGRAMS = dnsimp bug_1323 bug_58_double bug_79_double_complex bug_142 bug_142_gen
if ICB
check_PROGRAMS += icb_arpack_c
check_PROGRAMS += icb_arpack_cpp
check_PROGRAMS += bug_1315_single
check_PROGRAMS += bug_1315_double
endif
dist_check_DATA = testA.mtx
TESTS = $(check_PROGRAMS)

dnsimp_SOURCES = dnsimp.f mmio.f debug.h

bug_1323_SOURCES = bug_1323.f

bug_58_double_SOURCES = bug_58_double.f

bug_79_double_complex_SOURCES = bug_79_double_complex.f

bug_142_SOURCES = bug_142.f

bug_142_gen_SOURCES = bug_142_gen.f

if ICB
icb_arpack_c_SOURCES = icb_arpack_c.c
icb_arpack_c_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
icb_arpack_c_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir) -I$(top_srcdir)/ICB

icb_arpack_cpp_SOURCES = icb_arpack_cpp.cpp
icb_arpack_cpp_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
icb_arpack_cpp_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir) -I$(top_srcdir)/ICB

bug_1315_single_SOURCES = bug_1315_single.c
bug_1315_double_SOURCES = bug_1315_double.c
bug_1315_single_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
bug_1315_double_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
bug_1315_single_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir) -I$(top_srcdir)/ICB
bug_1315_double_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_builddir) -I$(top_srcdir)/ICB
endif
