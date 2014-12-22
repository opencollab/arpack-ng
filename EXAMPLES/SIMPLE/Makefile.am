LDADD = $(top_builddir)/libarpack.la $(BLAS_LIBS) $(LAPACK_LIBS)

SIMPLE = sssimp dssimp snsimp dnsimp cnsimp znsimp

check_PROGRAMS = $(SIMPLE)
TESTS = $(check_PROGRAMS)

sssimp_SOURCES = sssimp.f
dssimp_SOURCES = dssimp.f
snsimp_SOURCES = snsimp.f
dnsimp_SOURCES = dnsimp.f
cnsimp_SOURCES = cnsimp.f
znsimp_SOURCES = znsimp.f
