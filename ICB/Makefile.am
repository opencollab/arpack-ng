arpack_includedir = $(includedir)/arpack

if ICB
AM_DISTCHECK_CONFIGURE_FLAGS = --enable-icb
arpack_include_HEADERS = debug_c.h debug_c.hpp
arpack_include_HEADERS += stat_c.h stat_c.hpp
arpack_include_HEADERS += arpack.h arpack.hpp
if MPI
arpack_include_HEADERS += parpack.h parpack.hpp
endif
endif

noinst_LTLIBRARIES = libdbgicb.la libstaicb.la
libdbgicb_la_SOURCES = debug_icb.F90
libdbgicb_la_FCFLAGS = -I$(top_builddir)
libstaicb_la_SOURCES = stat_icb.F90
libstaicb_la_FCFLAGS = -I$(top_builddir)
