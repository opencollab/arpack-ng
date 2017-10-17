LDADD=$(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

check_PROGRAMS = dnsimp bug_1315_single bug_1315_double bug_1323 bug_58_double
if ICB
check_PROGRAMS += iso_c_binding_c
check_PROGRAMS += iso_c_binding_cpp
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

if ICB
iso_c_binding_c_SOURCES = iso_c_binding_c.c
iso_c_binding_c_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
iso_c_binding_c_CFLAGS = -I$(top_builddir)

iso_c_binding_cpp_SOURCES = iso_c_binding_cpp.cpp
iso_c_binding_cpp_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
iso_c_binding_cpp_CPPFLAGS = -I$(top_builddir)
endif
