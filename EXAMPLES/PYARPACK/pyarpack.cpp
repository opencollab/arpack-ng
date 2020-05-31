#include <Python.h>  // PyErr_SetString.

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <complex>
#include <pyarpackDrtSolver.hpp>
#include <pyarpackItrSolver.hpp>
#include <sstream>  // ostringstream.
#include <string>
#include <vector>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

template <typename RC, typename FD, typename EM, typename SLV>
void exportArpackSparseItr(bp::scope& pySlv, std::string const& dtype) {
  // Created nested namespace in module.

  pySlv.attr(dtype.c_str()) =
      bp::class_<pyarpackSparseItrSolver<RC, FD, EM, SLV>>(
          dtype.c_str(),
          "arpack data type (must be consistent with numpy dtype)")
          .def("solve", &pyarpackSparseItrSolver<RC, FD, EM, SLV>::solve,
               (bp::arg("A"), bp::arg("B") = bp::tuple()),
               "solve standard or generalised eigen problem where A and B must "
               "be sparse and provided in coo format: (dimension, row-indice "
               "array, column-indice array, matrice-value array) tuple")
          .def("checkEigVec",
               &pyarpackSparseItrSolver<RC, FD, EM, SLV>::checkEigVec,
               (bp::arg("A"), bp::arg("B") = bp::tuple(),
                bp::arg("diffTol") = 1.e-3),
               "check eigen vectors accuracy where A and B must be sparse and "
               "provided in coo format: (dimension, row-indice array, "
               "column-indice array, matrice-value array) tuple")
              ARPACKSOLVERMEMBER(pyarpackSparseItrSolver)
          .def_readwrite(
              "slvTol", &pyarpackSparseItrSolver<RC, FD, EM, SLV>::slvTol,
              "tolerance of the iterative mode solver - default: 1.e-6")
          .def_readwrite("slvMaxIt",
                         &pyarpackSparseItrSolver<RC, FD, EM, SLV>::slvMaxIt,
                         "maximum number of iterations of the iterative mode "
                         "solver - default: 100")
          .def_readwrite(
              "slvILUDropTol",
              &pyarpackSparseItrSolver<RC, FD, EM, SLV>::slvILUDropTol,
              "drop tolerance of the ILU preconditioner (if any) of the "
              "iterative mode solver - default: 1")
          .def_readwrite(
              "slvILUFillFactor",
              &pyarpackSparseItrSolver<RC, FD, EM, SLV>::slvILUFillFactor,
              "fill factor of the ILU preconditioner (if any) of the iterative "
              "mode solver - default: 2");
};

template <typename RC, typename FD, typename EM, typename SLV>
void exportArpackSparseDrt(bp::scope& pySlv, std::string const& dtype) {
  // Created nested namespace in module.

  pySlv.attr(dtype.c_str()) =
      bp::class_<pyarpackSparseDrtSolver<RC, FD, EM, SLV>>(
          dtype.c_str(),
          "arpack data type (must be consistent with numpy dtype)")
          .def("solve", &pyarpackSparseDrtSolver<RC, FD, EM, SLV>::solve,
               (bp::arg("A"), bp::arg("B") = bp::tuple()),
               "solve standard or generalised eigen problem where A and B must "
               "be sparse and provided in coo format: (dimension, row-indice "
               "array, column-indice array, matrice-value array) tuple")
          .def("checkEigVec",
               &pyarpackSparseDrtSolver<RC, FD, EM, SLV>::checkEigVec,
               (bp::arg("A"), bp::arg("B") = bp::tuple(),
                bp::arg("diffTol") = 1.e-3),
               "check eigen vectors accuracy where A and B must be sparse and "
               "provided in coo format: (dimension, row-indice array, "
               "column-indice array, matrice-value array) tuple")
              ARPACKSOLVERMEMBER(pyarpackSparseDrtSolver)
          .def_readwrite(
              "slvPvtThd", &pyarpackSparseDrtSolver<RC, FD, EM, SLV>::slvPvtThd,
              "pivoting tolerance of the direct mode solver - default: 1.e-6")
          .def_readwrite("slvOffset",
                         &pyarpackSparseDrtSolver<RC, FD, EM, SLV>::slvOffset,
                         "cholesky offset (LLT, LDLT) of the direct mode "
                         "solver - default: 0.")
          .def_readwrite("slvScale",
                         &pyarpackSparseDrtSolver<RC, FD, EM, SLV>::slvScale,
                         "cholesky scale (LLT, LDLT) of the direct mode solver "
                         "- default: 1.");
};

template <typename RC, typename FD, typename EM, typename SLV>
void exportArpackDenseDrt(bp::scope& pySlv, std::string const& dtype) {
  // Created nested namespace in module.

  pySlv.attr(dtype.c_str()) =
      bp::class_<pyarpackDenseDrtSolver<RC, FD, EM, SLV>>(
          dtype.c_str(),
          "arpack data type (must be consistent with numpy dtype)")
          .def("solve", &pyarpackDenseDrtSolver<RC, FD, EM, SLV>::solve,
               (bp::arg("A"), bp::arg("B") = bp::tuple()),
               "solve standard or generalised eigen problem where A and B must "
               "be dense and provided in raw format: (n-squared matrice-value "
               "array, row or column ordered boolean)")
          .def("checkEigVec",
               &pyarpackDenseDrtSolver<RC, FD, EM, SLV>::checkEigVec,
               (bp::arg("A"), bp::arg("B") = bp::tuple(),
                bp::arg("diffTol") = 1.e-3),
               "check eigen vectors accuracy where A and B must be dense and "
               "provided in raw format: (n-squared matrice-value array, row or "
               "column ordered boolean)")
              ARPACKSOLVERMEMBER(pyarpackDenseDrtSolver)
          .def_readwrite(
              "slvPvtThd", &pyarpackDenseDrtSolver<RC, FD, EM, SLV>::slvPvtThd,
              "pivoting tolerance of the direct mode solver - default: 1.e-6")
          .def_readwrite("slvOffset",
                         &pyarpackDenseDrtSolver<RC, FD, EM, SLV>::slvOffset,
                         "cholesky offset (LLT, LDLT) of the direct mode "
                         "solver - default: 0.")
          .def_readwrite("slvScale",
                         &pyarpackDenseDrtSolver<RC, FD, EM, SLV>::slvScale,
                         "cholesky scale (LLT, LDLT) of the direct mode solver "
                         "- default: 1.");
};

class sparseBiCGDiag {};
class sparseBiCGILU {};
class sparseCGDiag {};
class sparseCGILU {};
class sparseLLT {};
class sparseLDLT {};
class sparseLU {};
class sparseQR {};

class denseLLT {};
class denseLDLT {};
class denseLURR {};
class denseQRRR {};
class denseLUPP {};
class denseQRPP {};

std::complex<double> EigVecZGetItem(
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& M, int idx) {
  if (idx < 0 || idx >= M.size()) {
    pyarpackThrowError("index out of range");
    return std::complex<double>();
  }
  return M[idx];
};

std::string EigVecZToString(EigVecZ const& vec) {
  std::ostringstream s;
  s << vec;
  return s.str();
};

BOOST_PYTHON_MODULE(pyarpack) {
  // Initialize.

  bn::initialize();

  bp::class_<std::vector<std::complex<double>>>("StdVecZ").def(
      bp::vector_indexing_suite<std::vector<std::complex<double>>>());

  bp::class_<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>>("EigVecZ")
      .def("__getitem__", &EigVecZGetItem)
      .def("__str__", &EigVecZToString);

  bp::class_<std::vector<EigVecZ>>("StdVecEVZ")
      .def("__iter__", bp::iterator<std::vector<EigVecZ>>())
      .def(bp::vector_indexing_suite<std::vector<EigVecZ>>());

  // Documentation of the python module.

  std::ostringstream doc;
  doc << "You can use sparse or dense matrices, and, play with iterative or direct mode solvers (CG, LU, ...):" << std::endl;
  doc << "1. choose arpack solver with a given mode solver" << std::endl;
  doc << "   1.1. if you need to handle sparse matrices" << std::endl;
  doc << "        >> from pyarpack import sparseBiCG as pyarpackSlv" << std::endl;
  doc << "   1.2. if you need to handle dense matrices" << std::endl;
  doc << "        >> from pyarpack import denseBiCG as pyarpackSlv" << std::endl;
  doc << "2. choose arpack data type (float, double, ...)" << std::endl;
  doc << "   >> arpackSlv = pyarpackSlv.double()" << std::endl;
  doc << "3. solve the eigen problem" << std::endl;
  doc << "   >> arpackSlv.solve(A, B)" << std::endl;
  doc << "4. get eigen values and vectors" << std::endl;
  doc << "   >> print(arpackSlv.vec)" << std::endl;
  doc << "   >> print(arpackSlv.val)" << std::endl;
  doc << std::endl;
  doc << "Notes:" << std::endl;
  doc << "1. arpack data type (float, double, ...) must be consistent with A/B numpy dtypes (float32, float64, ...)." << std::endl;
  doc << "   at python side, the data MUST be casted in the EXACT expected type (int32, int64, float, double, ...)." << std::endl;
  doc << "   otherwise, C++ may not get the data the way it expects them: C++ will not know how to read python data." << std::endl;
  doc << "   if you are not sure how data have been passed from python to C++, set arpackSlv.debug = 1 and check out debug traces." << std::endl;
  doc << "   in other words, pyarpack users MUST :" << std::endl;
  doc << "   1.1. create numpy arrays specifying explicitly the type:" << std::endl;
  doc << "        >> Aij = np.array([], dtype='complex128')" << std::endl;
  doc << "   1.2. filling numpy arrays casting value on append:" << std::endl;
  doc << "        >> Aij = np.append(Aij, np.complex128(np.complex( 200.,  200.))) # Casting value on append is MANDATORY or C++ won't get the expected type." << std::endl;
  doc << "   1.3. calling the solver flavor which is consistent with the numpy array data type:" << std::endl;
  doc << "        >> arpackSlv = pyarpackSlv.complexDouble() # Caution: complexDouble <=> np.array(..., dtype='complex128')" << std::endl;
  doc << "   note: NO data type check can be done at C++ side, the pyarpack user MUST insure data consistency." << std::endl;
  doc << "2. sparse matrices must be provided in coo format (n, i, j, Mij), that is, as a tuple where:" << std::endl;
  doc << "   2.1. n is an integer." << std::endl;
  doc << "   2.2. i, j, Mij are 1 x nnz numpy arrays." << std::endl;
  doc << "3. dense  matrices must be provided in raw format (Mij, rowOrdered), that is, as a tuple where:" << std::endl;
  doc << "   3.1. Mij is an n x n numpy array." << std::endl;
  doc << "   3.2. rowOrdered is a boolean (column ordered if False)." << std::endl;
  doc << "4. arpack mode solver are provided by eigen:" << std::endl;
  doc << "   4.1. when solver is iterative, A and B can be sparse only." << std::endl;
  doc << "   4.2. when solver is direct,    A and B can be sparse or dense." << std::endl;
  bp::scope().attr("__doc__") = doc.str().c_str();

  // Specify that this module is actually a package.

  bp::object package = bp::scope();
  package.attr("__path__") = "pyarpack";

  // Create  python module.

  std::string module = "pyarpack";
  bp::object pyModule(
      bp::handle<>(bp::borrowed(PyImport_AddModule(module.c_str()))));

  // Create modules.

  {
    std::string slv = "sparseBiCGDiag";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): BiCG with diagonal (Jacobi) "
        "preconditioner";
    bp::scope pySlvBiCGDiag =
        bp::class_<sparseBiCGDiag>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseItr<float, float, EigSMxS, EigSBiCGS>(pySlvBiCGDiag,
                                                            "float");
    exportArpackSparseItr<double, double, EigSMxD, EigSBiCGD>(pySlvBiCGDiag,
                                                              "double");
    exportArpackSparseItr<std::complex<float>, float, EigSMxC, EigSBiCGC>(
        pySlvBiCGDiag, "complexFloat");
    exportArpackSparseItr<std::complex<double>, double, EigSMxZ, EigSBiCGZ>(
        pySlvBiCGDiag, "complexDouble");
  }

  {
    std::string slv = "sparseBiCGILU";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): BiCG with ILU preconditioner";
    bp::scope pySlvBiCGILU =
        bp::class_<sparseBiCGILU>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseItr<float, float, EigSMxS, EigSBiCGILUS>(pySlvBiCGILU,
                                                               "float");
    exportArpackSparseItr<double, double, EigSMxD, EigSBiCGILUD>(pySlvBiCGILU,
                                                                 "double");
    exportArpackSparseItr<std::complex<float>, float, EigSMxC, EigSBiCGILUC>(
        pySlvBiCGILU, "complexFloat");
    exportArpackSparseItr<std::complex<double>, double, EigSMxZ, EigSBiCGILUZ>(
        pySlvBiCGILU, "complexDouble");
  }

  {
    std::string slv = "sparseCGDiag";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): CG with diagonal (Jacobi) "
        "preconditioner";
    bp::scope pySlvCGDiag =
        bp::class_<sparseCGDiag>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseItr<float, float, EigSMxS, EigSCGS>(pySlvCGDiag, "float");
    exportArpackSparseItr<double, double, EigSMxD, EigSCGD>(pySlvCGDiag,
                                                            "double");
    exportArpackSparseItr<std::complex<float>, float, EigSMxC, EigSCGC>(
        pySlvCGDiag, "complexFloat");
    exportArpackSparseItr<std::complex<double>, double, EigSMxZ, EigSCGZ>(
        pySlvCGDiag, "complexDouble");
  }

  {
    std::string slv = "sparseCGILU";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): CG with ILU preconditioner";
    bp::scope pySlvCGILU =
        bp::class_<sparseCGILU>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseItr<float, float, EigSMxS, EigSCGILUS>(pySlvCGILU,
                                                             "float");
    exportArpackSparseItr<double, double, EigSMxD, EigSCGILUD>(pySlvCGILU,
                                                               "double");
    exportArpackSparseItr<std::complex<float>, float, EigSMxC, EigSCGILUC>(
        pySlvCGILU, "complexFloat");
    exportArpackSparseItr<std::complex<double>, double, EigSMxZ, EigSCGILUZ>(
        pySlvCGILU, "complexDouble");
  }

  {
    std::string slv = "sparseLLT";
    std::string slvHelp = "arpack internal mode solver (mode > 1): LLT";
    bp::scope pySlvLLT = bp::class_<sparseLLT>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseDrt<float, float, EigSMxS, EigSLLTS>(pySlvLLT, "float");
    exportArpackSparseDrt<double, double, EigSMxD, EigSLLTD>(pySlvLLT,
                                                             "double");
    exportArpackSparseDrt<std::complex<float>, float, EigSMxC, EigSLLTC>(
        pySlvLLT, "complexFloat");
    exportArpackSparseDrt<std::complex<double>, double, EigSMxZ, EigSLLTZ>(
        pySlvLLT, "complexDouble");
  }

  {
    std::string slv = "sparseLDLT";
    std::string slvHelp = "arpack internal mode solver (mode > 1): LDLT";
    bp::scope pySlvLDLT = bp::class_<sparseLDLT>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseDrt<float, float, EigSMxS, EigSLDLTS>(pySlvLDLT, "float");
    exportArpackSparseDrt<double, double, EigSMxD, EigSLDLTD>(pySlvLDLT,
                                                              "double");
    exportArpackSparseDrt<std::complex<float>, float, EigSMxC, EigSLDLTC>(
        pySlvLDLT, "complexFloat");
    exportArpackSparseDrt<std::complex<double>, double, EigSMxZ, EigSLDLTZ>(
        pySlvLDLT, "complexDouble");
  }

  {
    std::string slv = "sparseLU";
    std::string slvHelp = "arpack internal mode solver (mode > 1): LU";
    bp::scope pySlvLU = bp::class_<sparseLU>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseDrt<float, float, EigSMxS, EigSLUS>(pySlvLU, "float");
    exportArpackSparseDrt<double, double, EigSMxD, EigSLUD>(pySlvLU, "double");
    exportArpackSparseDrt<std::complex<float>, float, EigSMxC, EigSLUC>(
        pySlvLU, "complexFloat");
    exportArpackSparseDrt<std::complex<double>, double, EigSMxZ, EigSLUZ>(
        pySlvLU, "complexDouble");
  }

  {
    std::string slv = "sparseQR";
    std::string slvHelp = "arpack internal mode solver (mode > 1): QR";
    bp::scope pySlvQR = bp::class_<sparseQR>(slv.c_str(), slvHelp.c_str());
    exportArpackSparseDrt<float, float, EigSMxS, EigSQRS>(pySlvQR, "float");
    exportArpackSparseDrt<double, double, EigSMxD, EigSQRD>(pySlvQR, "double");
    exportArpackSparseDrt<std::complex<float>, float, EigSMxC, EigSQRC>(
        pySlvQR, "complexFloat");
    exportArpackSparseDrt<std::complex<double>, double, EigSMxZ, EigSQRZ>(
        pySlvQR, "complexDouble");
  }

  {
    std::string slv = "denseLLT";
    std::string slvHelp = "arpack internal mode solver (mode > 1): LLT";
    bp::scope pySlvLLT = bp::class_<denseLLT>(slv.c_str(), slvHelp.c_str());
    exportArpackDenseDrt<float, float, EigDMxS, EigDLLTS>(pySlvLLT, "float");
    exportArpackDenseDrt<double, double, EigDMxD, EigDLLTD>(pySlvLLT, "double");
    exportArpackDenseDrt<std::complex<float>, float, EigDMxC, EigDLLTC>(
        pySlvLLT, "complexFloat");
    exportArpackDenseDrt<std::complex<double>, double, EigDMxZ, EigDLLTZ>(
        pySlvLLT, "complexDouble");
  }

  {
    std::string slv = "denseLDLT";
    std::string slvHelp = "arpack internal mode solver (mode > 1): LDLT";
    bp::scope pySlvLDLT = bp::class_<denseLDLT>(slv.c_str(), slvHelp.c_str());
    exportArpackDenseDrt<float, float, EigDMxS, EigDLDLTS>(pySlvLDLT, "float");
    exportArpackDenseDrt<double, double, EigDMxD, EigDLDLTD>(pySlvLDLT,
                                                             "double");
    exportArpackDenseDrt<std::complex<float>, float, EigDMxC, EigDLDLTC>(
        pySlvLDLT, "complexFloat");
    exportArpackDenseDrt<std::complex<double>, double, EigDMxZ, EigDLDLTZ>(
        pySlvLDLT, "complexDouble");
  }

  {
    std::string slv = "denseLURR";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): LU Rank Revealing (slower, "
        "more stable)";
    bp::scope pySlvLURR = bp::class_<denseLURR>(slv.c_str(), slvHelp.c_str());
    exportArpackDenseDrt<float, float, EigDMxS, EigDFLUS>(pySlvLURR, "float");
    exportArpackDenseDrt<double, double, EigDMxD, EigDFLUD>(pySlvLURR,
                                                            "double");
    exportArpackDenseDrt<std::complex<float>, float, EigDMxC, EigDFLUC>(
        pySlvLURR, "complexFloat");
    exportArpackDenseDrt<std::complex<double>, double, EigDMxZ, EigDFLUZ>(
        pySlvLURR, "complexDouble");
  }

  {
    std::string slv = "denseQRRR";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): QR Rank Revealing (slower, "
        "more stable)";
    bp::scope pySlvQRRR = bp::class_<denseQRRR>(slv.c_str(), slvHelp.c_str());
    exportArpackDenseDrt<float, float, EigDMxS, EigDFQRS>(pySlvQRRR, "float");
    exportArpackDenseDrt<double, double, EigDMxD, EigDFQRD>(pySlvQRRR,
                                                            "double");
    exportArpackDenseDrt<std::complex<float>, float, EigDMxC, EigDFQRC>(
        pySlvQRRR, "complexFloat");
    exportArpackDenseDrt<std::complex<double>, double, EigDMxZ, EigDFQRZ>(
        pySlvQRRR, "complexDouble");
  }

  {
    std::string slv = "denseLUPP";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): LU Partial Pivoting (faster, "
        "less stable)";
    bp::scope pySlvLUPP = bp::class_<denseLUPP>(slv.c_str(), slvHelp.c_str());
    exportArpackDenseDrt<float, float, EigDMxS, EigDPLUS>(pySlvLUPP, "float");
    exportArpackDenseDrt<double, double, EigDMxD, EigDPLUD>(pySlvLUPP,
                                                            "double");
    exportArpackDenseDrt<std::complex<float>, float, EigDMxC, EigDPLUC>(
        pySlvLUPP, "complexFloat");
    exportArpackDenseDrt<std::complex<double>, double, EigDMxZ, EigDPLUZ>(
        pySlvLUPP, "complexDouble");
  }

  {
    std::string slv = "denseQRPP";
    std::string slvHelp =
        "arpack internal mode solver (mode > 1): QR Partial Pivoting (faster, "
        "less stable)";
    bp::scope pySlvQPPR = bp::class_<denseQRPP>(slv.c_str(), slvHelp.c_str());
    exportArpackDenseDrt<float, float, EigDMxS, EigDPQRS>(pySlvQPPR, "float");
    exportArpackDenseDrt<double, double, EigDMxD, EigDPQRD>(pySlvQPPR,
                                                            "double");
    exportArpackDenseDrt<std::complex<float>, float, EigDMxC, EigDPQRC>(
        pySlvQPPR, "complexFloat");
    exportArpackDenseDrt<std::complex<double>, double, EigDMxZ, EigDPQRZ>(
        pySlvQPPR, "complexDouble");
  }
}

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=2 ts=2 et smartindent :*/
