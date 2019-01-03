LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

EXTRA_DIST = README

check_PROGRAMS = arpackmm
dist_check_SCRIPTS = arpackmm.sh
dist_check_DATA = An.mtx As.mtx B.mtx
TESTS = arpackmm.sh

arpackmm_SOURCES = arpackmm.cpp
arpackmm_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
arpackmm_CPPFLAGS = -I$(top_srcdir)/ICB $(EIGEN3_CFLAGS)

CLEANFILES = resid.out v.out
