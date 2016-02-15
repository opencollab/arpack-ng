LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

SIMPLE = sssimp dssimp snsimp dnsimp cnsimp znsimp

check_PROGRAMS = $(SIMPLE)
TESTS = $(check_PROGRAMS)

EXTRA_DIST = README debug.h

sssimp_SOURCES = sssimp.f
dssimp_SOURCES = dssimp.f
snsimp_SOURCES = snsimp.f
dnsimp_SOURCES = dnsimp.f
cnsimp_SOURCES = cnsimp.f
znsimp_SOURCES = znsimp.f
