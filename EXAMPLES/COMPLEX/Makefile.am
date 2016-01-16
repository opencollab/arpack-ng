LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

CNDRV = cndrv1 cndrv2 cndrv3 cndrv4
ZNDRV = zndrv1 zndrv2 zndrv3 zndrv4
COMPLEX = $(CNDRV) $(ZNDRV)

check_PROGRAMS = $(COMPLEX)
TESTS = $(check_PROGRAMS)

EXTRA_DIST = README

# Complex problem using single complex
cndrv1_SOURCES = cndrv1.f
cndrv2_SOURCES = cndrv2.f
cndrv3_SOURCES = cndrv3.f
cndrv4_SOURCES = cndrv4.f

# Complex problem using double complex
zndrv1_SOURCES = zndrv1.f
zndrv2_SOURCES = zndrv2.f
zndrv3_SOURCES = zndrv3.f
zndrv4_SOURCES = zndrv4.f
