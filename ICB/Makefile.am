AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)

pkgincludedir = $(includedir)/arpack

pkginclude_HEADERS = debug_c.h debug_c.hpp
pkginclude_HEADERS += stat_c.h stat_c.hpp
pkginclude_HEADERS += iseed_c.h iseed_c.hpp
pkginclude_HEADERS += arpack.h arpack.hpp
if MPI
pkginclude_HEADERS += parpack.h parpack.hpp
endif

noinst_LTLIBRARIES = libdbgicb.la libstaicb.la libisdicb.la
libdbgicb_la_SOURCES = debug_icb.F90
libstaicb_la_SOURCES = stat_icb.F90
libisdicb_la_SOURCES = iseed_icb.F90
