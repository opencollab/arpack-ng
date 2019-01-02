pkgincludedir = $(includedir)/arpack

if ICB
AM_DISTCHECK_CONFIGURE_FLAGS = --enable-icb
pkginclude_HEADERS = debug_c.h debug_c.hpp
pkginclude_HEADERS += stat_c.h stat_c.hpp
pkginclude_HEADERS += arpack.h arpack.hpp
if MPI
pkginclude_HEADERS += parpack.h parpack.hpp
endif
endif

noinst_LTLIBRARIES = libdbgicb.la libstaicb.la
libdbgicb_la_SOURCES = debug_icb.F90
libdbgicb_la_FCFLAGS = -I$(top_builddir)
libstaicb_la_SOURCES = stat_icb.F90
libstaicb_la_FCFLAGS = -I$(top_builddir)
