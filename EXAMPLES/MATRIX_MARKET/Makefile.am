LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX)$(ITF64SUFFIX).la
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/ICB $(EIGEN3_CFLAGS)

EXTRA_DIST = README

check_PROGRAMS = arpackmm
dist_check_SCRIPTS = arpackmm.sh issue401.sh issue215.sh
dist_check_DATA = An.mtx As.mtx Az.mtx B.mtx Bz.mtx issue401.mtx issue215.mtx
TESTS = arpackmm.sh issue401.sh issue215.sh

arpackmm_SOURCES = arpackmm.cpp

pkgincludedir = $(includedir)/arpack@ITF64SUFFIX@
if EIGEN
pkginclude_HEADERS = arpackSolver.hpp
endif

CLEANFILES = \
	arpackmm.run.log \
	arpackmm.resid.out \
	arpackmm.v.out \
	arpackSolver.resid.out \
	arpackSolver.v.out
