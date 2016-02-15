LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

SSDRV = ssdrv1 ssdrv2 ssdrv3 ssdrv4 ssdrv5 ssdrv6
DSDRV = dsdrv1 dsdrv2 dsdrv3 dsdrv4 dsdrv5 dsdrv6
SYM = $(SSDRV) $(DSDRV)

check_PROGRAMS = $(SYM)
TESTS = $(check_PROGRAMS)

EXTRA_DIST = README

# simple symmetric problem using single precision
ssdrv1_SOURCES = ssdrv1.f
ssdrv2_SOURCES = ssdrv2.f
ssdrv3_SOURCES = ssdrv3.f
ssdrv4_SOURCES = ssdrv4.f
ssdrv5_SOURCES = ssdrv5.f
ssdrv6_SOURCES = ssdrv6.f

# simple symmetric problem using double precision
dsdrv1_SOURCES = dsdrv1.f
dsdrv2_SOURCES = dsdrv2.f
dsdrv3_SOURCES = dsdrv3.f
dsdrv4_SOURCES = dsdrv4.f
dsdrv5_SOURCES = dsdrv5.f
dsdrv6_SOURCES = dsdrv6.f
