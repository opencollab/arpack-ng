LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/ICB $(EIGEN3_CFLAGS)

EXTRA_DIST = README

check_PROGRAMS = arpackmm
dist_check_SCRIPTS = arpackmm.sh
dist_check_DATA = An.mtx As.mtx Az.mtx B.mtx Bz.mtx
SH_LOG_COMPILER = bash
TESTS = arpackmm.sh

arpackmm_SOURCES = arpackmm.cpp

CLEANFILES = resid.out v.out
