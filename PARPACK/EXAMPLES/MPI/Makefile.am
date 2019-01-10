F77 = $(MPIF77)
LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

SNDRV = psndrv1 psndrv3
DNDRV = pdndrv1 pdndrv3
SSDRV = pssdrv1
DSDRV = pdsdrv1
CNDRV = pcndrv1
ZNDRV = pzndrv1

check_PROGRAMS = $(SNDRV) $(DNDRV) $(SSDRV) $(DSDRV) $(CNDRV) $(ZNDRV)

# Run MPI tests with "mpirun -n 2"
LOG_COMPILER = mpirun
LOG_FLAGS = -n 2

TESTS = $(check_PROGRAMS)

EXTRA_DIST = debug.h stat.h README.CALLING-PARPACK-FROM-C-OR-CPP

# Simple nonsymmetric problem using single precision
psndrv1_SOURCES = psndrv1.f
psndrv3_SOURCES = psndrv3.f

# Simple nonsymmetric problem using double precision
pdndrv1_SOURCES = pdndrv1.f
pdndrv3_SOURCES = pdndrv3.f

# Simple symmetric problem using single precision
pssdrv1_SOURCES = pssdrv1.f

# Simple symmetric problem using double precision
pdsdrv1_SOURCES = pdsdrv1.f

# Complex problem using single complex
pcndrv1_SOURCES = pcndrv1.f

# Complex problem using double complex
pzndrv1_SOURCES = pzndrv1.f
