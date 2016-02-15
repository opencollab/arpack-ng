LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

SNDRV = sndrv1 sndrv2 sndrv3 sndrv4 sndrv5 sndrv6
DNDRV = dndrv1 dndrv2 dndrv3 dndrv4 dndrv5 dndrv6
NONSYM = $(SNDRV) $(DNDRV)

check_PROGRAMS = $(NONSYM)
TESTS = $(check_PROGRAMS)

EXTRA_DIST = README

# Simple nonsymmetric problem using single precision
sndrv1_SOURCES = sndrv1.f
sndrv2_SOURCES = sndrv2.f
sndrv3_SOURCES = sndrv3.f
sndrv4_SOURCES = sndrv4.f
sndrv5_SOURCES = sndrv5.f
sndrv6_SOURCES = sndrv6.f

# Simple nonsymmetric problem using double precision
dndrv1_SOURCES = dndrv1.f
dndrv2_SOURCES = dndrv2.f
dndrv3_SOURCES = dndrv3.f
dndrv4_SOURCES = dndrv4.f
dndrv5_SOURCES = dndrv5.f
dndrv6_SOURCES = dndrv6.f
