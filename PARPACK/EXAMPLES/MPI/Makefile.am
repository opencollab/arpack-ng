F77 = $(MPIF77)

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
psndrv1_SOURCES = psndrv1.F90
psndrv1_CPPFLAGS = -I$(top_srcdir)
psndrv1_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)
psndrv3_SOURCES = psndrv3.F90
psndrv3_CPPFLAGS = -I$(top_srcdir)
psndrv3_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)

# Simple nonsymmetric problem using double precision
pdndrv1_SOURCES = pdndrv1.F90
pdndrv1_CPPFLAGS = -I$(top_srcdir)
pdndrv1_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)
pdndrv3_SOURCES = pdndrv3.F90
pdndrv3_CPPFLAGS = -I$(top_srcdir)
pdndrv3_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)

# Simple symmetric problem using single precision
pssdrv1_SOURCES = pssdrv1.F90
pssdrv1_CPPFLAGS = -I$(top_srcdir)
pssdrv1_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)

# Simple symmetric problem using double precision
pdsdrv1_SOURCES = pdsdrv1.F90
pdsdrv1_CPPFLAGS = -I$(top_srcdir)
pdsdrv1_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)

# Complex problem using single complex
pcndrv1_SOURCES = pcndrv1.F90
pcndrv1_CPPFLAGS = -I$(top_srcdir)
pcndrv1_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)

# Complex problem using double complex
pzndrv1_SOURCES = pzndrv1.F90
pzndrv1_CPPFLAGS = -I$(top_srcdir)
pzndrv1_LDADD = $(top_builddir)/PARPACK/SRC/MPI/libparpack@LIBSUFFIX@@ITF64SUFFIX@.la $(LAPACK_LIBS) $(BLAS_LIBS)
