LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/ICB $(EIGEN3_CFLAGS)

EXTRA_DIST = README

lib_LTLIBRARIES = libpyarpack@LIBSUFFIX@.la
libpyarpack@LIBSUFFIX@_la_SOURCES = pyarpack.cpp
libpyarpack@LIBSUFFIX@_la_CPPFLAGS = $(BOOST_CPPFLAGS) $(PYTHON_CPPFLAGS) -I$(top_builddir)/EXAMPLES/MATRIX_MARKET -I$(top_builddir)/ICB -I$(top_builddir) $(EIGEN3_CFLAGS)
libpyarpack@LIBSUFFIX@_la_LDFLAGS = $(BOOST_LDFLAGS) $(PYTHON_LDFLAGS)
libpyarpack@LIBSUFFIX@_la_LIBADD = $(BOOST_NUMPY_LIB) $(BOOST_PYTHON_LIB) $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS)

SH_LOG_COMPILER = $(PYTHON)
TESTS  = pyarpackDenseLDLT.py
TESTS += pyarpackDenseLLT.py
TESTS += pyarpackDenseLUPP.py
TESTS += pyarpackDenseLURR.py
TESTS += pyarpackDenseQRPP.py
TESTS += pyarpackDenseQRRR.py
TESTS += pyarpackSparseBiCGDiag.py
TESTS += pyarpackSparseBiCGILU.py
TESTS += pyarpackSparseCGDiag.py
TESTS += pyarpackSparseCGILU.py
TESTS += pyarpackSparseLDLT.py
TESTS += pyarpackSparseLLT.py
TESTS += pyarpackSparseLU.py
TESTS += pyarpackSparseQR.py
