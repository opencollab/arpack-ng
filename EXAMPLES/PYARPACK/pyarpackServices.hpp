#ifndef __PYARPACKSERVICES_HPP__
#define __PYARPACKSERVICES_HPP__

#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include <cmath> // sqrt.

#include <Eigen/Sparse>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

#define ARPACKSOLVERMEMBER(pyarpackSolver)                                                          \
    .def_readwrite("symPb",           &pyarpackSolver<RC, FD, EM, SLV>::symPb,                      \
                   "symmetric problem - default: true")                                             \
    .def_readwrite("nbEV",            &pyarpackSolver<RC, FD, EM, SLV>::nbEV,                       \
                   "number of eigen vectors to find - default: 1")                                  \
    .def_readwrite("nbCV",            &pyarpackSolver<RC, FD, EM, SLV>::nbCV,                       \
                   "number of dimensions of the workspace - default: 3")                            \
    .def_readwrite("tol",             &pyarpackSolver<RC, FD, EM, SLV>::tol,                        \
                   "tolerance - default: 1.e-6")                                                    \
    .def_readwrite("sigmaReal",       &pyarpackSolver<RC, FD, EM, SLV>::sigmaReal,                  \
                   "shift over real axis - default: 0.")                                            \
    .def_readwrite("sigmaImag",       &pyarpackSolver<RC, FD, EM, SLV>::sigmaImag,                  \
                   "shift over imaginary axis - default: 0.")                                       \
    .def_readwrite("dumpToFile",      &pyarpackSolver<RC, FD, EM, SLV>::dumpToFile,                 \
                   "dump eigen vectors to arpackSolver.*.out files - default: false")               \
    .def_readwrite("restartFromFile", &pyarpackSolver<RC, FD, EM, SLV>::restartFromFile,            \
                   "restart from eigen vectors found in arpackSolver.*.out files - default: false") \
    .def_readwrite("mag",             &pyarpackSolver<RC, FD, EM, SLV>::mag,                        \
                   "magnitude - default: LM")                                                       \
    .def_readwrite("maxIt",           &pyarpackSolver<RC, FD, EM, SLV>::maxIt,                      \
                   "maximum number of arpack iterations - default: 100")                            \
    .def_readwrite("schur",           &pyarpackSolver<RC, FD, EM, SLV>::schur,                      \
                   "compute schur vectors - default: false")                                        \
    .def_readwrite("verbose",         &pyarpackSolver<RC, FD, EM, SLV>::verbose,                    \
                   "verbosity level - default: 0")                                                  \
    .def_readonly ("stdPb",           &pyarpackSolver<RC, FD, EM, SLV>::stdPb,                      \
                   "standard or generalised problem - default: true")                               \
    .def_readonly ("val",             &pyarpackSolver<RC, FD, EM, SLV>::val,                        \
                   "eigen values found")                                                            \
    .def_readonly ("vec",             &pyarpackSolver<RC, FD, EM, SLV>::vec,                        \
                   "eigen vectors found")                                                           \
    .def_readonly ("mode",            &pyarpackSolver<RC, FD, EM, SLV>::mode,                       \
                   "selected arpack mode (according to input options: std/gen, shift, ...)")        \
    .def_readonly ("nbIt",            &pyarpackSolver<RC, FD, EM, SLV>::nbIt,                       \
                   "number of arpack iterations")                                                   \
    .def_readonly ("imsTime",         &pyarpackSolver<RC, FD, EM, SLV>::imsTime,                    \
                   "time spent to initialize the mode solver if needed")                            \
    .def_readonly ("rciTime",         &pyarpackSolver<RC, FD, EM, SLV>::rciTime,                    \
                   "time spent in Reverse Communication Interface")                                 \
    .def_readwrite("debug",           &pyarpackSolver<RC, FD, EM, SLV>::debug,                      \
                   "debug traces (up to 3) - default: 0")                                           \

#define ARPACKSOLVERDEBUGSTAT()                                                             \
if (debug > 3) debug = 3;                                                                   \
debug_c(6, -6, debug, debug, debug, debug, debug, debug, debug, debug, debug, debug, debug, \
        debug, debug, debug, debug, debug, debug, debug, debug, debug, debug, debug);       \
stat_c(nopx, nbx, nrorth, nitref, nrstrt, tsaupd, tsaup2,                                   \
       tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd, tnaup2,                              \
       tnaitr, tneigt, tngets, tnapps, tnconv, tcaupd, tcaup2,                              \
       tcaitr, tceigt, tcgets, tcapps, tcconv, tmvopx, tmvbx,                               \
       tgetv0, titref, trvec);                                                              \

void pyarpackThrowError(std::string const & msg) {
  std::string const info = "Error: " + msg;
  std::cerr << info << std::endl;
  PyErr_SetString(PyExc_IndexError, info.c_str());
  bp::throw_error_already_set();
};

template<typename RC, typename EM>
class pyarpackServices {
  // Public methods.

  public:

    static int buildSparseMatrice(bp::tuple const & T, Eigen::SparseMatrix<RC> & M,
                                  a_int const & debug, std::string const & msg) {
      // Get boost data as C++ data.

      if (bp::len(T) != 4) {pyarpackThrowError(msg + " must be a 3-tuple"); return 1;}
      bp::extract<int>           nExt(T[0]);
      bp::extract<bn::ndarray>   iExt(T[1]);
      bp::extract<bn::ndarray>   jExt(T[2]);
      bp::extract<bn::ndarray> mijExt(T[3]);
      if (!  nExt.check()) {pyarpackThrowError(msg + "[0] must be an integer" ); return 1;}
      if (!  iExt.check()) {pyarpackThrowError(msg + "[1] must be numpy.array"); return 1;}
      if (!  jExt.check()) {pyarpackThrowError(msg + "[2] must be numpy.array"); return 1;}
      if (!mijExt.check()) {pyarpackThrowError(msg + "[3] must be numpy.array"); return 1;}
      bn::ndarray   iArray =   iExt();
      bn::ndarray   jArray =   jExt();
      bn::ndarray mijArray = mijExt();
      if (iArray.get_dtype()   != bn::dtype::get_builtin<a_int>()) {pyarpackThrowError(msg + "[1] type is not consistent"); return 1;}
      if (jArray.get_dtype()   != bn::dtype::get_builtin<a_int>()) {pyarpackThrowError(msg + "[2] type is not consistent"); return 1;}
      if (mijArray.get_dtype() != bn::dtype::get_builtin<RC>()   ) {pyarpackThrowError(msg + "[3] type is not consistent with arpack type"); return 1;}

      a_int iSz = iArray.shape(0);
      a_int * iPtr = reinterpret_cast<a_int*>(iArray.get_data());
      a_int jSz = jArray.shape(0);
      a_int * jPtr = reinterpret_cast<a_int*>(jArray.get_data());
      a_int mSz = mijArray.shape(0);
      RC * mPtr = reinterpret_cast<RC*>(mijArray.get_data());

      if (iSz != jSz) {pyarpackThrowError(msg + "[1] and " + msg + "[2] must have same lenght"); return 1;}
      if (iSz != mSz) {pyarpackThrowError(msg + "[1] and " + msg + "[3] must have same lenght"); return 1;}

      // Debug on demand: casting value on numpy.append is MANDATORY or C++ won't get the expected type..

      for (auto k = 0; debug && k < mSz; k++) {
        std::cout << "pyarpackServices::buildSparseMatrice - " << msg << "[" << iPtr[k] << ", " << jPtr[k] << "] = " << mPtr[k] << std::endl;
      };

      // Build sparse matrice.

      a_uint n = nExt();
      a_uint iMin = n+1, jMin = n+1;
      for (auto k = 0; k < mSz; k++) {
        if (iPtr[k] < iMin) iMin = iPtr[k];
        if (jPtr[k] < jMin) jMin = jPtr[k];
      };
      if (iMin != 0 && iMin != 1) {pyarpackThrowError(msg + ": smallest row indice must be 0 or 1"); return 1;}
      if (jMin != 0 && jMin != 1) {pyarpackThrowError(msg + ": smallest column indice must be 0 or 1"); return 1;}
      a_int iBased = 0, jBased = 0;
      if (iMin == 1) iBased = 1;
      if (jMin == 1) jBased = 1;

      M = Eigen::SparseMatrix<RC>(n, n); // Set matrice dimensions.
      std::vector<Eigen::Triplet<RC>> triplets;
      a_uint nnz = mSz;
      triplets.reserve(nnz);
      for (auto k = 0; k < nnz; k++) triplets.emplace_back(iPtr[k] - iBased, jPtr[k] - jBased, mPtr[k]);
      M.setFromTriplets(triplets.begin(), triplets.end()); // Set all (i, j, Mij).

      return 0;
    };

    static int buildDenseMatrice(bp::tuple const & T, Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> & M,
                                 a_int const & debug, std::string const & msg) {
      // Get boost data as C++ data.

      if (bp::len(T) != 2) {pyarpackThrowError(msg + " must be a 2-tuple"); return 1;}
      bp::extract<bn::ndarray> mijExt(T[0]);
      bp::extract<bool>          oExt(T[1]);
      if (!mijExt.check()) {pyarpackThrowError(msg + " must be numpy.array"); return 1;}
      if (  !oExt.check()) {pyarpackThrowError(msg + " must be a boolean"); return 1;}
      bn::ndarray   mijArray = mijExt();
      bool        rowOrdered = oExt();
      if (mijArray.get_dtype() != bn::dtype::get_builtin<RC>()) {pyarpackThrowError(msg + " type is not consistent with arpack type"); return 1;}

      a_int mSz = mijArray.shape(0);
      RC * mPtr = reinterpret_cast<RC*>(mijArray.get_data());

      a_uint n = std::sqrt(mSz);
      if (n*n != mSz) {pyarpackThrowError(msg + " must be a squared matrice"); return 1;}

      // Debug on demand: casting value on numpy.append is MANDATORY or C++ won't get the expected type..

      for (auto k = 0; debug && k < mSz; k++) {
        std::cout << "pyarpackServices::buildDenseMatrice - " << msg << "[" << k << "] = " << mPtr[k] << std::endl;
      };

      // Build dense matrice.

      M = Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic>(n, n); // Set matrice dimensions.
      M.setZero(n, n); // Avoid spurious/random values which may break solves (LU, QR, ...).
      if (rowOrdered) {
        for (size_t k = 0; k < n; k++) {
          for (size_t l = 0; l < n; l++) M(k, l) = mPtr[l+k*n];
        }
      }
      else {
        for (size_t l = 0; l < n; l++) {
          for (size_t k = 0; k < n; k++) M(k, l) = mPtr[k+l*n];
        }
      }

      return 0;
    };
};

#endif

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=2 ts=2 et smartindent :*/
