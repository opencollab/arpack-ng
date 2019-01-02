LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

EXTRA_DIST = README

check_PROGRAMS = arpackmm
TESTS = arpackmm.sh

arpackmm_SOURCES = arpackmm.cpp
arpackmm_LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
arpackmm_CPPFLAGS = -I$(top_srcdir)/ICB $(EIGEN3_CFLAGS)
