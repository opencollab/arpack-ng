// This code sample is meant for convenience (not performance):
//   - test/run arpack (eigen values / vectors, timing).
//   - play with modes: shift, invert, shift + invert.
//   - use with user matrices (matrix market format).

#include <iostream>
#include <string>
#include <sstream> // stringstream.
#include <fstream> // [io]fstream.
#include <vector>
#include <complex>
#include <algorithm> // max_element.
#include <chrono>
#include <limits> // epsilon.
#include <cmath> // fabs.
#include <iomanip> // setw.
#include <cassert> // assert.
#include <memory> // unique_ptr.
#include "arpack.h"
#include "debug_c.hpp"
#include "stat_c.hpp"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Cholesky>

using namespace std;

// Vector related types.

typedef Eigen::Matrix<         float,  Eigen::Dynamic, 1> EigVecS; // Real.
typedef Eigen::Matrix<        double,  Eigen::Dynamic, 1> EigVecD; // Real.
typedef Eigen::Matrix<complex< float>, Eigen::Dynamic, 1> EigVecC; // Complex.
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> EigVecZ; // Complex.

typedef Eigen::Map<EigVecS> EigMpVS; // Real.
typedef Eigen::Map<EigVecD> EigMpVD; // Real.
typedef Eigen::Map<EigVecC> EigMpVC; // Complex.
typedef Eigen::Map<EigVecZ> EigMpVZ; // Complex.

// Sparse matrix related types.

typedef Eigen::SparseMatrix<         float>  EigSMxS; // Real.
typedef Eigen::Triplet     <         float>  EigSMTS; // Real.
typedef Eigen::SparseMatrix<        double>  EigSMxD; // Real.
typedef Eigen::Triplet     <        double>  EigSMTD; // Real.
typedef Eigen::SparseMatrix<complex< float>> EigSMxC; // Complex.
typedef Eigen::Triplet     <complex< float>> EigSMTC; // Complex.
typedef Eigen::SparseMatrix<complex<double>> EigSMxZ; // Complex.
typedef Eigen::Triplet     <complex<double>> EigSMTZ; // Complex.

typedef Eigen::BiCGSTAB         <EigSMxS> EigSBiCGS; // Real.
typedef Eigen::BiCGSTAB         <EigSMxD> EigSBiCGD; // Real.
typedef Eigen::BiCGSTAB         <EigSMxC> EigSBiCGC; // Complex.
typedef Eigen::BiCGSTAB         <EigSMxZ> EigSBiCGZ; // Complex.

typedef Eigen::ConjugateGradient<EigSMxS> EigSCGS;   // Real.
typedef Eigen::ConjugateGradient<EigSMxD> EigSCGD;   // Real.
typedef Eigen::ConjugateGradient<EigSMxC> EigSCGC;   // Complex.
typedef Eigen::ConjugateGradient<EigSMxZ> EigSCGZ;   // Complex.

typedef Eigen::IncompleteLUT<         float>                                  EigILUS;     // Real.
typedef Eigen::IncompleteLUT<        double>                                  EigILUD;     // Real.
typedef Eigen::IncompleteLUT<complex< float>>                                 EigILUC;     // Complex.
typedef Eigen::IncompleteLUT<complex<double>>                                 EigILUZ;     // Complex.

typedef Eigen::BiCGSTAB         <EigSMxS,                            EigILUS> EigSBiCGILUS; // Real.
typedef Eigen::BiCGSTAB         <EigSMxD,                            EigILUD> EigSBiCGILUD; // Real.
typedef Eigen::BiCGSTAB         <EigSMxC,                            EigILUC> EigSBiCGILUC; // Complex.
typedef Eigen::BiCGSTAB         <EigSMxZ,                            EigILUZ> EigSBiCGILUZ; // Complex.

typedef Eigen::ConjugateGradient<EigSMxS, Eigen::Lower|Eigen::Upper, EigILUS> EigSCGILUS;   // Real.
typedef Eigen::ConjugateGradient<EigSMxD, Eigen::Lower|Eigen::Upper, EigILUD> EigSCGILUD;   // Real.
typedef Eigen::ConjugateGradient<EigSMxC, Eigen::Lower|Eigen::Upper, EigILUC> EigSCGILUC;   // Complex.
typedef Eigen::ConjugateGradient<EigSMxZ, Eigen::Lower|Eigen::Upper, EigILUZ> EigSCGILUZ;   // Complex.

typedef Eigen::SimplicialLLT <EigSMxS, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLLTS;  // Real.
typedef Eigen::SimplicialLLT <EigSMxD, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLLTD;  // Real.
typedef Eigen::SimplicialLLT <EigSMxC, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLLTC;  // Complex.
typedef Eigen::SimplicialLLT <EigSMxZ, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLLTZ;  // Complex.

typedef Eigen::SimplicialLDLT<EigSMxS, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLDLTS; // Real.
typedef Eigen::SimplicialLDLT<EigSMxD, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLDLTD; // Real.
typedef Eigen::SimplicialLDLT<EigSMxC, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLDLTC; // Complex.
typedef Eigen::SimplicialLDLT<EigSMxZ, Eigen::Lower, Eigen::COLAMDOrdering<int>> EigSLDLTZ; // Complex.

typedef Eigen::SparseLU<EigSMxS, Eigen::COLAMDOrdering<int>> EigSLUS;  // Real.
typedef Eigen::SparseLU<EigSMxD, Eigen::COLAMDOrdering<int>> EigSLUD;  // Real.
typedef Eigen::SparseLU<EigSMxC, Eigen::COLAMDOrdering<int>> EigSLUC;  // Complex.
typedef Eigen::SparseLU<EigSMxZ, Eigen::COLAMDOrdering<int>> EigSLUZ;  // Complex.

typedef Eigen::SparseQR<EigSMxS, Eigen::COLAMDOrdering<int>> EigSQRS;  // Real.
typedef Eigen::SparseQR<EigSMxD, Eigen::COLAMDOrdering<int>> EigSQRD;  // Real.
typedef Eigen::SparseQR<EigSMxC, Eigen::COLAMDOrdering<int>> EigSQRC;  // Complex.
typedef Eigen::SparseQR<EigSMxZ, Eigen::COLAMDOrdering<int>> EigSQRZ;  // Complex.

// Dense matrix related types.

typedef Eigen::Matrix<         float,  Eigen::Dynamic, Eigen::Dynamic> EigDMxS; // Real.
typedef Eigen::Matrix<        double,  Eigen::Dynamic, Eigen::Dynamic> EigDMxD; // Real.
typedef Eigen::Matrix<complex< float>, Eigen::Dynamic, Eigen::Dynamic> EigDMxC; // Complex.
typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> EigDMxZ; // Complex.

typedef Eigen::LLT<EigDMxS> EigDLLTS; // Real.
typedef Eigen::LLT<EigDMxD> EigDLLTD; // Real.
typedef Eigen::LLT<EigDMxC> EigDLLTC; // Complex.
typedef Eigen::LLT<EigDMxZ> EigDLLTZ; // Complex.

typedef Eigen::LDLT<EigDMxS> EigDLDLTS; // Real.
typedef Eigen::LDLT<EigDMxD> EigDLDLTD; // Real.
typedef Eigen::LDLT<EigDMxC> EigDLDLTC; // Complex.
typedef Eigen::LDLT<EigDMxZ> EigDLDLTZ; // Complex.

typedef Eigen::FullPivLU<EigDMxS> EigDFLUS; // Real.
typedef Eigen::FullPivLU<EigDMxD> EigDFLUD; // Real.
typedef Eigen::FullPivLU<EigDMxC> EigDFLUC; // Complex.
typedef Eigen::FullPivLU<EigDMxZ> EigDFLUZ; // Complex.

typedef Eigen::FullPivHouseholderQR<EigDMxS> EigDFQRS; // Real.
typedef Eigen::FullPivHouseholderQR<EigDMxD> EigDFQRD; // Real.
typedef Eigen::FullPivHouseholderQR<EigDMxC> EigDFQRC; // Complex.
typedef Eigen::FullPivHouseholderQR<EigDMxZ> EigDFQRZ; // Complex.

typedef Eigen::PartialPivLU<EigDMxS> EigDPLUS; // Real.
typedef Eigen::PartialPivLU<EigDMxD> EigDPLUD; // Real.
typedef Eigen::PartialPivLU<EigDMxC> EigDPLUC; // Complex.
typedef Eigen::PartialPivLU<EigDMxZ> EigDPLUZ; // Complex.

typedef Eigen::HouseholderQR<EigDMxS> EigDPQRS; // Real.
typedef Eigen::HouseholderQR<EigDMxD> EigDPQRD; // Real.
typedef Eigen::HouseholderQR<EigDMxC> EigDPQRC; // Complex.
typedef Eigen::HouseholderQR<EigDMxZ> EigDPQRZ; // Complex.

// Class, functions, program.

class options {
  public:
    options() {
      fileA = "A.mtx";
      fileB = "N.A."; // Not available.
      dense = false;
      denseRR = true;
      nbEV = 1;
      nbCV = 2*nbEV + 1;
      stdPb = true; // Standard or generalized (= not standard).
      symPb = true;
      cpxPb = false;
      simplePrec = false; // Double precision.
      mag = string("LM"); // Large magnitude.
      shiftReal = false; shiftImag = false;
      sigmaReal = 0.; sigmaImag = 0.; // Eigen value translation: look for lambda+sigma instead of lambda.
      invert = false; // Eigen value invertion: look for 1./lambda instead of lambda.
      tol = 1.e-06;
      maxIt = 100;
      schur = false; // Compute Ritz vectors.
      slv = "BiCG";
      slvItrTol = nullptr;
      slvItrMaxIt = nullptr;
      slvItrPC = "diag";
      check = true;
      verbose = 0;
      debug = 0;
      restart = false;
    };

    int readCmdLine(int argc, char ** argv) {
      // Check for command line independent parameters.

      for (int a = 1; argv && a < argc; a++) {
        string clo = argv[a]; // Command line option.
        if (clo == "--help") return usage(0);
        if (clo == "--A") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          fileA = argv[a];
        }
        if (clo == "--dense") {
          dense = true;
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          string rr(argv[a]);
          if (rr != "true" && rr != "false") {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
          denseRR = (rr == "true") ? true : false;
        }
        if (clo == "--nbEV") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream nEV(argv[a]);
          nEV >> nbEV; if (!nEV) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
          nbCV = 2*nbEV + 1;
        }
        if (clo == "--genPb") {
          stdPb = false;
          fileB = "B.mtx";
        }
        if (clo == "--nonSymPb") symPb = false;
        if (clo == "--cpxPb") {
          symPb = false;
          cpxPb = true;
        }
        if (clo == "--simplePrec") simplePrec = true;
        if (clo == "--mag") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          mag = argv[a]; // small mag (likely poor perf) <=> large mag + invert (likely good perf).
          bool ok = (mag == "LM" || mag == "SM" || mag == "LR" || mag == "SR" || mag == "LI" || mag == "SI") ? true : false;
          if (!ok) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--shiftReal") {
          shiftReal = true;
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream s(argv[a]);
          s >> sigmaReal; if (!s) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--shiftImag") {
          shiftImag = true;
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream s(argv[a]);
          s >> sigmaImag; if (!s) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--invert") invert = true;
        if (clo == "--tol") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream t(argv[a]);
          t >> tol; if (!t) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--maxIt") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream mi(argv[a]);
          mi >> maxIt; if (!mi) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--slv") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          slv = argv[a];
        }
        if (clo == "--slvItrTol") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream t(argv[a]);
          double tol = 0.;
          t >> tol; if (!t) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
          slvItrTol = unique_ptr<double>(new double);
          if (slvItrTol) *slvItrTol = tol;
        }
        if (clo == "--slvItrMaxIt") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream mi(argv[a]);
          int maxIt = 0;
          mi >> maxIt; if (!mi) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
          slvItrMaxIt = unique_ptr<int>(new int);
          if (slvItrMaxIt) *slvItrMaxIt = maxIt;
        }
        if (clo == "--slvItrPC") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream pc(argv[a]);
          pc >> slvItrPC; if (!pc) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--noCheck") check = false;
        if (clo == "--verbose") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream vb(argv[a]);
          vb >> verbose; if (!vb) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--debug") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream dbg(argv[a]);
          dbg >> debug; if (!dbg) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
          if (debug > 3) debug = 3;
          debug_c(6, -6, debug, debug, debug, debug, debug, debug, debug, debug, debug, debug, debug,
                  debug, debug, debug, debug, debug, debug, debug, debug, debug, debug, debug);
        }
        if (clo == "--restart") restart = true;
      }

      // Check for command line dependent parameters.

      for (int a = 1; argv && a < argc; a++) {
        string clo = argv[a]; // Command line option.
        if (clo == "--nbCV") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream nCV(argv[a]);
          nCV >> nbCV; if (!nCV) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--B") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          fileB = argv[a];
        }
      }

      return 0;
    };

    int usage(int rc = 1) {
      cout << "Usage: running arpack with matrix market files to check for eigen values/vectors." << endl;
      cout << endl;
      cout << "  --A F:            file name of matrix A such that A X = lambda X. (standard)" << endl;
      cout << "                    the file F must be compliant with the matrix market format." << endl;
      cout << "                    default: A.mtx" << endl;
      cout << "  --B F:            file name of matrix B such that A X = lambda B X. (generalized)" << endl;
      cout << "                    the file F must be compliant with the matrix market format." << endl;
      cout << "                    default: N.A. for standard problem, or, B.mtx for generalized problem" << endl;
      cout << "  --dense RR:       consider A and B as dense matrices." << endl;
      cout << "                    if RR = true,  use more-stable-but-slow versions of LU / QR (rank revealing)." << endl;
      cout << "                    if RR = false, use less-stable-but-fast versions of LU / QR (depends on condition number)." << endl;
      cout << "                    Notes:" << endl;
      cout << "                      - only direct solvers are available when using dense matrices." << endl;
      cout << "                      - RR does not impact the use of LLT and LDLT." << endl;
      cout << "                      - thresholds only make sense for rank-revealing decompositions." << endl;
      cout << "                    default: consider A and B as sparse matrices" << endl;
      cout << "  --nbEV:           number of eigen values/vectors to compute." << endl;
      cout << "                    default: 1" << endl;
      cout << "  --nbCV:           number of columns of the matrix V." << endl;
      cout << "                    default: 2*nbEV+1" << endl;
      cout << "  --genPb:          generalized problem." << endl;
      cout << "                    default: standard problem" << endl;
      cout << "  --nonSymPb:       non symmetric problem (<=> use dn[ae]upd)." << endl;
      cout << "                    default: symmetric problem (<=> use ds[ae]upd)" << endl;
      cout << "  --cpxPb:          complex (non symmetric) problem (<=> use zn[ae]upd)." << endl;
      cout << "                    default: false (<=> use d*[ae]upd)" << endl;
      cout << "  --simplePrec:     use simple precision (less accurate, but, half memory footprint)." << endl;
      cout << "                    default: false (<=> use double precision: use [dz]*upd)" << endl;
      cout << "  --mag M:          set magnitude of eigen values to look for (LM, SM, LR, SR, LI, SI)." << endl;
      cout << "                    default: large magnitude (LM)" << endl;
      cout << "  --shiftReal S:    real shift where sigma = S (look for lambda+S instead of lambda)." << endl;
      cout << "                    default: no shift, S = 0." << endl;
      cout << "  --shiftImag S:    imaginary shift where sigma = S (look for lambda+S instead of lambda)." << endl;
      cout << "                    default: no shift, S = 0." << endl;
      cout << "  --invert:         invert mode (look for 1./lambda instead of lambda)." << endl;
      cout << "                    default: no invert" << endl;
      cout << "  --tol T:          tolerance T." << endl;
      cout << "                    default: 1.e-06" << endl;
      cout << "  --maxIt M:        maximum iterations M." << endl;
      cout << "                    default: 100" << endl;
      cout << "  --schur:          compute Schur vectors." << endl;
      cout << "                      the Schur decomposition is such that A = Q^H x T x Q where:" << endl;
      cout << "                        - the H superscript refers to the Hermitian transpose: Q^H = (Q^t)^*." << endl;
      cout << "                        - Q is unitary: Q is such that Q^H x Q = I." << endl;
      cout << "                        - T is an upper-triangular matrix whose diagonal elements are the eigenvalues of A." << endl;
      cout << "                      every square matrix has a Schur decomposition: columns of Q are the Schur vectors." << endl;
      cout << "                      for a general matrix A, there is no relation between Schur vectors of A and eigenvectors of A." << endl;
      cout << "                      if q_j is the i-th Schur vector, then A x q_j is a linear combination of q_1, ..., q_j." << endl;
      cout << "                      Schur vectors q_1, q_2, ..., q_j span an invariant subspace of A." << endl;
      cout << "                      the Schur vectors and eigenvectors of A are the same if A is a normal matrix." << endl;
      cout << "                    default: compute Ritz vectors (approximations of eigen vectors)" << endl;
      cout << "  --slv S:          solver (needed if arpack mode > 1)." << endl;
      cout << "                      BiCG:     iterative method, any matrices" << endl;
      cout << "                      CG:       iterative method, sym matrices only" << endl;
      cout << "                      LU#P:     direct method, any matrices (pivoting needed)" << endl;
      cout << "                        P:        pivoting threshold (not used if not RR)" << endl;
      cout << "                      QR#P:     direct method, any matrices (pivoting needed)" << endl;
      cout << "                        P:        pivoting threshold (not used if not RR)" << endl;
      cout << "                      LLT#O#S:  direct method, SPD matrices only (pivoting not needed)" << endl;
      cout << "                        O:        shift offset (not used if --dense)" << endl;
      cout << "                        S:        shift scale (not used if --dense)" << endl;
      cout << "                      LDLT#O#S: direct method, symmetric positive semi-definite matrices only (pivoting not needed)" << endl;
      cout << "                        O:        shift offset (not used if --dense)" << endl;
      cout << "                        S:        shift scale (not used if --dense)" << endl;
      cout << "                    default: BiCG" << endl;
      cout << "  --slvItrTol T:    solver tolerance T (for iterative solvers)." << endl;
      cout << "                    default: eigen default value" << endl;
      cout << "  --slvItrMaxIt M:  solver maximum iterations M (for iterative solvers)." << endl;
      cout << "                    default: eigen default value" << endl;
      cout << "  --slvItrPC PC:    solver preconditioner (for iterative solvers)." << endl;
      cout << "                      PC preconditioner:" << endl;
      cout << "                        diag:    eigen diagonal preconditioner (Jacobi)." << endl;
      cout << "                        ILU#D#F: eigen ILU preconditioner." << endl;
      cout << "                          D:       drop tolerance." << endl;
      cout << "                          F:       fill factor." << endl;
      cout << "                    default: diagonal preconditioner (Jacobi)" << endl;
      cout << "  --noCheck:        check arpack eigen values/vectors." << endl;
      cout << "                    check will fail if Schur vectors are computed and A is NOT a normal matrix." << endl;
      cout << "                    default: check" << endl;
      cout << "  --verbose V:      verbosity level (up to 3)." << endl;
      cout << "                    default: 0" << endl;
      cout << "  --debug D:        debug level (up to 3)." << endl;
      cout << "                    default: 0" << endl;
      cout << "  --restart:        restart from previous run (which had produced arpackmm.*.out)." << endl;
      cout << "                    restart from eigen basis approximation computed during a previous run." << endl;
      cout << "                    default: false" << endl;
      if (rc == 0) exit(0);
      return rc;
    };

    friend ostream & operator<< (ostream & ostr, options const & opt);

    string fileA;
    string fileB;
    bool dense;
    bool denseRR;
    a_int nbEV;
    a_int nbCV;
    bool stdPb; // Standard or generalized (= not standard).
    bool symPb;
    bool cpxPb;
    bool simplePrec;
    string mag; // Magnitude <=> "which" arpack parameter.
    bool shiftReal, shiftImag;
    double sigmaReal, sigmaImag; // Eigen value translation: look for lambda+sigma instead of lambda.
    bool invert; // Eigen value invertion: look for 1./lambda instead of lambda.
    double tol;
    int maxIt;
    bool schur;
    string slv;
    unique_ptr<double> slvItrTol;
    unique_ptr<int> slvItrMaxIt;
    string slvItrPC;
    bool check;
    int verbose;
    int debug;
    bool restart;
};

ostream & operator<< (ostream & ostr, options const & opt) {
  ostr << "OPT: A " << opt.fileA << ", B " << opt.fileB;
  if      (opt.dense &&  opt.denseRR) ostr << ", dense yes (RR true)";
  else if (opt.dense && !opt.denseRR) ostr << ", dense yes (RR false)";
  else                                ostr << ", dense no";
  ostr << ", nbEV " << opt.nbEV << ", nbCV " << opt.nbCV << ", stdPb " << (opt.stdPb ? "yes" : "no");
  ostr << ", symPb " << (opt.symPb ? "yes" : "no") <<  ", cpxPb " << (opt.cpxPb ? "yes" : "no");
  ostr <<  ", simplePrec " << (opt.simplePrec ? "yes" : "no") << ", mag " << opt.mag << endl;
  ostr << "OPT: shiftReal " << (opt.shiftReal ? "yes" : "no") << ", sigmaReal " << opt.sigmaReal;
  ostr << ", shiftImag " << (opt.shiftImag ? "yes" : "no") << ", sigmaImag " << opt.sigmaImag;
  ostr << ", invert " << (opt.invert ? "yes" : "no") << ", tol " << opt.tol << ", maxIt " << opt.maxIt;
  ostr << ", " << (opt.schur ? "Schur" : "Ritz") << " vectors" << endl;
  ostr << "OPT: slv " << opt.slv << ", slvItrPC " << opt.slvItrPC;
  if (opt.slvItrTol)    ostr << ", slvItrTol "    << *opt.slvItrTol;
  if (opt.slvItrMaxIt)  ostr << ", slvItrMaxIt "  << *opt.slvItrMaxIt;
  ostr << ", check " << (opt.check ? "yes" : "no") << ", verbose " << opt.verbose << ", debug " << opt.debug;
  ostr << ", restart " << (opt.restart ? "yes" : "no") << endl;
  return ostr;
}

void makeZero(         float  & zero) {zero = 0.f;}

void makeZero(        double  & zero) {zero = 0.;}

void makeZero(complex< float> & zero) {zero = complex<double>(0.f, 0.f);}

void makeZero(complex<double> & zero) {zero = complex<double>(0., 0.);}

template<typename RC>
int readMatrixMarket(string const & fileName,
                     a_uint & n, a_uint & m, vector<a_uint> & i, vector<a_uint> & j, vector<RC> & Mij) {
  ifstream inp(fileName);
  if (!inp) {cerr << "Error: can not open " << fileName << endl; return 1;}

  // Read matrix from file.

  a_uint l = 0, nnz = 0;
  do {
    // Skip comments.

    string inpLine; getline(inp, inpLine); l++;
    while (isspace(*inpLine.begin())) inpLine.erase(inpLine.begin()); // Suppress leading white spaces.
    if (inpLine.length() == 0) continue; // Empty line.
    if (inpLine[0] == '%') continue; // Comments skipped, begin reading.

    // Read matrix market file.

    stringstream inpSS(inpLine);
    if (n == 0 && m == 0) { // Header.
      inpSS >> n >> m;
      if (!inpSS) {cerr << "Error: bad header (n, m)" << endl; return 1;}
      if (nnz == 0) {
        inpSS >> nnz;
        if (inpSS) { // OK, (optional) nnz has been provided.
          i.reserve(nnz);
          j.reserve(nnz);
          Mij.reserve(nnz);
        }
      }
    }
    else { // Body.
      a_uint k = 0, l = 0;
      RC zero; makeZero(zero);
      RC Mkl = zero;
      inpSS >> k >> l >> Mkl;
      if (!inpSS) {cerr << "Error: bad line (" << fileName << ", line " << l << ")" << endl; return 1;}
      i.push_back(k);
      j.push_back(l);
      Mij.push_back(Mkl);
    }
  }
  while (inp);

  // Handle 1-based -> 0-based.

  nnz = i.size(); // In case nnz was not provided.
  if (*max_element(begin(i), end(i)) == n || *max_element(begin(j), end(j)) == m) {
    for (size_t k = 0; k < nnz; k++) i[k] -= 1;
    for (size_t k = 0; k < nnz; k++) j[k] -= 1;
  }

  return 0;
}

template<typename RC, typename ET>
int createMatrix(string const & fileName, Eigen::SparseMatrix<RC> & M,
                 int const & verbose, string const & msg) {
  // Read matrix from file.

  a_uint n = 0, m = 0;
  vector<a_uint> i, j;
  vector<RC> Mij;
  int rc = readMatrixMarket<RC>(fileName, n, m, i, j, Mij);
  if (rc != 0) {cerr << "Error: read matrix market file KO" << endl; return rc;}

  // Create matrix from file.

  M = Eigen::SparseMatrix<RC>(n, m); // Set matrice dimensions.
  vector<ET> triplets;
  a_uint nnz = Mij.size();
  triplets.reserve(nnz);
  for (size_t k = 0; k < nnz; k++) triplets.emplace_back(i[k], j[k], Mij[k]);
  M.setFromTriplets(triplets.begin(), triplets.end()); // Set all (i, j, Mij).

  if (verbose == 3) {
    cout << endl << msg << endl;
    cout << endl << M << endl;
  }

  return 0;
}

template<typename RC, typename ET>
int createMatrix(string const & fileName, Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> & M,
                 int const & verbose, string const & msg) {
  // Read matrix from file.

  a_uint n = 0, m = 0;
  vector<a_uint> i, j;
  vector<RC> Mij;
  int rc = readMatrixMarket<RC>(fileName, n, m, i, j, Mij);
  if (rc != 0) {cerr << "Error: read matrix market file KO" << endl; return rc;}

  // Create matrix from file.

  M = Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic>(n, m); // Set matrice dimensions.
  M.setZero(n, m); // Avoid spurious/random values which may break solves (LU, QR, ...).
  a_uint nnz = Mij.size();
  for (size_t k = 0; k < nnz; k++) M(i[k], j[k]) = Mij[k];

  if (verbose == 3) {
    cout << endl << msg << endl;
    cout << endl << M << endl;
  }

  return 0;
}

class arpackEV { // Arpack eigen values / vectors.
  public:
    arpackEV() {mode = 1; nbIt = 0; imsTime = 0.; rciTime = 0.;};
    vector<complex<double>> val; // Eigen values.
    vector<EigVecZ> vec; // Eigen vectors.
    int mode;
    int nbIt;
    double imsTime; // Init mode solver.
    double rciTime;
};

void arpackAUPD(options const & opt,
                a_int * ido, char const * bMat, a_int nbDim, char const * which, float * resid, float * v,
                a_int ldv, a_int * iparam, a_int * ipntr, float * workd, float * workl, a_int lworkl, float * & rwork,
                a_int * info) {
  assert(rwork == NULL);
  if (opt.symPb) {
    ssaupd_c(ido, bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
  }
  else {
    snaupd_c(ido, bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
  }
}

void arpackAUPD(options const & opt,
                a_int * ido, char const * bMat, a_int nbDim, char const * which, double * resid, double * v,
                a_int ldv, a_int * iparam, a_int * ipntr, double * workd, double * workl, a_int lworkl, double * & rwork,
                a_int * info) {
  assert(rwork == NULL);
  if (opt.symPb) {
    dsaupd_c(ido, bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
  }
  else {
    dnaupd_c(ido, bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
  }
}

void arpackAUPD(options const & opt,
                a_int * ido, char const * bMat, a_int nbDim, char const * which, complex<float> * resid, complex<float> * v,
                a_int ldv, a_int * iparam, a_int * ipntr, complex<float> * workd, complex<float> * workl, a_int lworkl, float * & rwork,
                a_int * info) {
  if (!rwork) rwork = new float[opt.nbCV];
  cnaupd_c(ido, bMat, nbDim, which, opt.nbEV, opt.tol, reinterpret_cast<_Complex float*>(resid), opt.nbCV,
           reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr, reinterpret_cast<_Complex float*>(workd),
           reinterpret_cast<_Complex float*>(workl), lworkl, rwork, info);
}

void arpackAUPD(options const & opt,
                a_int * ido, char const * bMat, a_int nbDim, char const * which, complex<double> * resid, complex<double> * v,
                a_int ldv, a_int * iparam, a_int * ipntr, complex<double> * workd, complex<double> * workl, a_int lworkl, double * & rwork,
                a_int * info) {
  if (!rwork) rwork = new double[opt.nbCV];
  znaupd_c(ido, bMat, nbDim, which, opt.nbEV, opt.tol, reinterpret_cast<_Complex double*>(resid), opt.nbCV,
           reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr, reinterpret_cast<_Complex double*>(workd),
           reinterpret_cast<_Complex double*>(workl), lworkl, rwork, info);
}

template<typename RC, typename EV, typename EMV>
void arpackSpectrum(RC * d, RC * z, a_int nbDim, a_int * iparam, options const & opt, arpackEV & out) {
  // Arpack compute the whole spectrum.

  a_int nbConv = iparam[4];
  out.val.reserve(nbConv);
  for (a_int i = 0; d && i < nbConv; i++) {
    complex<double> lambda(d[i]);
    out.val.push_back(lambda);
    if (out.val.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).
  }

  out.vec.reserve(nbConv);
  for (a_int i = 0; z && i < nbConv; i++) {
    EV V = EMV(z + i*nbDim, nbDim);
    out.vec.push_back(V.template cast<complex<double>>());
    if (out.vec.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).
  }
}

template<typename RC, typename EV, typename EMV>
void arpackHalfSpectrum(RC * dr, RC * di, RC * z, a_int nbDim, a_int * iparam, options const & opt, arpackEV & out) {
  // Arpack compute only half of the spectrum.

  a_int nbConv = iparam[4];
  out.val.reserve(nbConv);
  for (a_int i = 0; dr && di && i <= nbConv/2; i++) { // Scan first half of the spectrum.
    // Get first half of the spectrum.

    complex<double> lambda(dr[i], di[i]);
    out.val.push_back(lambda);
    if (out.val.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).

    // Deduce second half of the spectrum.

    out.val.push_back(complex<double>(lambda.real(), -1.*lambda.imag()));
    if (out.val.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).
  }

  out.vec.reserve(nbConv);
  for (a_int i = 0; z && i <= nbConv/2; i++) { // Scan half spectrum.
    // Get first half of the spectrum.

    EV Vr = EMV(z + (2*i+0)*nbDim, nbDim); // Real part.
    EV Vi = EMV(z + (2*i+1)*nbDim, nbDim); // Imaginary part.
    complex<double> imag(0., 1.);
    EigVecZ V = Vr.template cast<complex<double>>() + imag * Vi.template cast<complex<double>>();
    out.vec.push_back(V);
    if (out.vec.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).

    // Deduce second half of the spectrum.

    V = Vr.template cast<complex<double>>() - imag * Vi.template cast<complex<double>>();
    out.vec.push_back(V);
    if (out.vec.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).
  }
}

int arpackEUPD(options const & opt, arpackEV & out,
               bool rvec, char const * howmny, a_int const * select, float * z,
               a_int ldz, char const * bMat, a_int nbDim, char const * which, float * resid, float * v,
               a_int ldv, a_int * iparam, a_int * ipntr, float * workd, float * workl, a_int lworkl, float * rwork,
               a_int & info) {
  assert(rwork == NULL);
  if (opt.symPb) {
    float * d = new float[opt.nbEV]; for (a_int k = 0; k < opt.nbEV; k++) d[k] = 0.;

    sseupd_c(rvec, howmny, select, d, z, ldz, opt.sigmaReal,
             bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
    if (info == -14) cerr << "Error: dseupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
    if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dseupd - KO with info " << info << endl; return 1;}

    arpackSpectrum<float, EigVecS, EigMpVS>(d, z, nbDim, iparam, opt, out);

    if (d) {delete [] d; d = NULL;}
  }
  else {
    float * dr = new float[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) dr[k] = 0.;
    float * di = new float[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) di[k] = 0.;
    float * workev = new float[3*opt.nbCV];

    sneupd_c(rvec, howmny, select, dr, di, z, ldz, opt.sigmaReal, opt.sigmaImag, workev,
             bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
    if (info == -14) cerr << "Error: dneupd - KO: [dz]naupd did not find any eigenvalues to sufficient accuracy" << endl;
    if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dneupd - KO with info " << info << endl; return 1;}

    arpackHalfSpectrum<float, EigVecS, EigMpVS>(dr, di, z, nbDim, iparam, opt, out);

    if (workev) {delete [] workev; workev = NULL;}
    if (dr)     {delete [] dr;     dr     = NULL;}
    if (di)     {delete [] di;     di     = NULL;}
  }

  return 0;
}

int arpackEUPD(options const & opt, arpackEV & out,
               bool rvec, char const * howmny, a_int const * select, double * z,
               a_int ldz, char const * bMat, a_int nbDim, char const * which, double * resid, double * v,
               a_int ldv, a_int * iparam, a_int * ipntr, double * workd, double * workl, a_int lworkl, double * rwork,
               a_int & info) {
  assert(rwork == NULL);
  if (opt.symPb) {
    double * d = new double[opt.nbEV]; for (a_int k = 0; k < opt.nbEV; k++) d[k] = 0.;

    dseupd_c(rvec, howmny, select, d, z, ldz, opt.sigmaReal,
             bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
    if (info == -14) cerr << "Error: dseupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
    if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dseupd - KO with info " << info << endl; return 1;}

    arpackSpectrum<double, EigVecD, EigMpVD>(d, z, nbDim, iparam, opt, out);

    if (d) {delete [] d; d = NULL;}
  }
  else {
    double * dr = new double[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) dr[k] = 0.;
    double * di = new double[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) di[k] = 0.;
    double * workev = new double[3*opt.nbCV];

    dneupd_c(rvec, howmny, select, dr, di, z, ldz, opt.sigmaReal, opt.sigmaImag, workev,
             bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
    if (info == -14) cerr << "Error: dneupd - KO: [dz]naupd did not find any eigenvalues to sufficient accuracy" << endl;
    if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dneupd - KO with info " << info << endl; return 1;}

    arpackHalfSpectrum<double, EigVecD, EigMpVD>(dr, di, z, nbDim, iparam, opt, out);

    if (workev) {delete [] workev; workev = NULL;}
    if (dr)     {delete [] dr;     dr     = NULL;}
    if (di)     {delete [] di;     di     = NULL;}
  }

  return 0;
}

int arpackEUPD(options const & opt, arpackEV & out,
               bool rvec, char const * howmny, a_int const * select, complex<float> * z,
               a_int ldz, char const * bMat, a_int nbDim, char const * which, complex<float> * resid, complex<float> * v,
               a_int ldv, a_int * iparam, a_int * ipntr, complex<float> * workd, complex<float> * workl, a_int lworkl, float * rwork,
               a_int & info) {
  complex<float> * d = new complex<float>[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) d[k] = complex<float>(0., 0.);
  complex<float> * workev = new complex<float>[2*opt.nbCV];

  complex<float> sigma = complex<float>((float) opt.sigmaReal, (float) opt.sigmaImag);
  cneupd_c(rvec, howmny, select, reinterpret_cast<_Complex float*>(d), reinterpret_cast<_Complex float*>(z), ldz,
           reinterpret_cast<_Complex float &>(sigma), reinterpret_cast<_Complex float*>(workev),
           bMat, nbDim, which, opt.nbEV, opt.tol, reinterpret_cast<_Complex float*>(resid), opt.nbCV,
           reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
           reinterpret_cast<_Complex float*>(workd), reinterpret_cast<_Complex float*>(workl), lworkl, rwork, &info);
  if (info == -14) cerr << "Error: zneupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
  if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: zneupd - KO with info " << info << endl; return 1;}

  arpackSpectrum<complex<float>, EigVecC, EigMpVC>(d, z, nbDim, iparam, opt, out);

  if (workev) {delete [] workev; workev = NULL;}
  if (d)      {delete [] d;           d = NULL;}

  return 0;
}

int arpackEUPD(options const & opt, arpackEV & out,
               bool rvec, char const * howmny, a_int const * select, complex<double> * z,
               a_int ldz, char const * bMat, a_int nbDim, char const * which, complex<double> * resid, complex<double> * v,
               a_int ldv, a_int * iparam, a_int * ipntr, complex<double> * workd, complex<double> * workl, a_int lworkl, double * rwork,
               a_int & info) {
  complex<double> * d = new complex<double>[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) d[k] = complex<double>(0., 0.);
  complex<double> * workev = new complex<double>[2*opt.nbCV];

  complex<double> sigma = complex<double>(opt.sigmaReal, opt.sigmaImag);
  zneupd_c(rvec, howmny, select, reinterpret_cast<_Complex double*>(d), reinterpret_cast<_Complex double*>(z), ldz,
           reinterpret_cast<_Complex double &>(sigma), reinterpret_cast<_Complex double*>(workev),
           bMat, nbDim, which, opt.nbEV, opt.tol, reinterpret_cast<_Complex double*>(resid), opt.nbCV,
           reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
           reinterpret_cast<_Complex double*>(workd), reinterpret_cast<_Complex double*>(workl), lworkl, rwork, &info);
  if (info == -14) cerr << "Error: zneupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
  if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: zneupd - KO with info " << info << endl; return 1;}

  arpackSpectrum<complex<double>, EigVecZ, EigMpVZ>(d, z, nbDim, iparam, opt, out);

  if (workev) {delete [] workev; workev = NULL;}
  if (d)      {delete [] d;           d = NULL;}

  return 0;
}

template<typename SLV, typename RC, typename FD, typename EM>
int arpackMode(options const & opt, int const mode,
               EM const & A, EM const & B, SLV & solver) {
  int rc = 1;

  if (mode == 1) {
    rc = 0;
  }
  else if (mode == 2 || mode == 3) {
    if (mode == 2) { // Invert mode.
      solver.compute(B);
    }
    else { // Shift invert mode.
      RC sigma; makeSigma(opt, sigma);
      auto S = A - sigma*B;
      solver.compute(S);
    }

    rc = 0;
  }
  else {cerr << "Error: arpack mode must be 1, 2 or 3 - KO" << endl; rc = 1;}

  return rc;
}

template<typename SLV, typename RC, typename FD>
int arpackMode(options const & opt, int const mode,
               Eigen::SparseMatrix<RC> const & A,
               Eigen::SparseMatrix<RC> const & B,
               SLV & solver) {
  return arpackMode<SLV, RC, FD, Eigen::SparseMatrix<RC>>(opt, mode, A, B, solver);
}

template<typename SLV, typename RC, typename FD>
int arpackMode(options const & opt, int const mode,
               Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> const & A,
               Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> & B,
               SLV & solver) {
  return arpackMode<SLV, RC, FD, Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic>>(opt, mode, A, B, solver);
}

template<typename RC>
void arpackRestart(options const & opt, string const & fileName,
                   a_int const & nbDim, RC * rv, bool readNbCV) {
  if (opt.restart) {
    ifstream ifs(fileName.c_str());
    if (ifs.is_open()) {
      a_int nbCV = 1; if (readNbCV) {ifs >> nbCV; if (opt.nbCV < nbCV) nbCV = opt.nbCV;}
      for (a_int n = 0; rv && n < nbDim*nbCV; n++) ifs >> rv[n];
      if (opt.verbose >= 1) {
        cout << endl;
        cout << fileName << ": restart OK" << endl;
        if (opt.verbose >= 2) {
          for (a_int n = 0; rv && n < nbDim; n++) cout << rv[n] << endl;
        }
        cout << endl;
      }
    }
  }
}

template<typename RC, typename FD,
         typename EM, typename EV,
         typename SLV>
int arpackSolve(options const & opt,
                EM const & A, EM const & B, SLV & solver, arpackEV & out) {
  // Arpack set up.

  // Note: some in/out parameters passed to d[sn][ae]upd are set to 0. before use.
  // d[sn][ae]upd uses dgetv0 to generate a random starting vector (when info is initialized to 0).
  // dgetv0 rely on resid/v: resid/v should be initialized to 0.0 to avoid "bad" starting random vectors.

  char const * which = opt.mag.c_str();
  a_int ido = 0; // First call to arpack.
  char const * iMat = "I";
  char const * gMat = "G";
  char const * bMat = (out.mode == 1) ? iMat : gMat;
  a_int nbDim = A.rows();
  RC zero; makeZero(zero);
  RC * resid = new RC[nbDim]; for (a_int n = 0; n < nbDim; n++) resid[n] = zero; // Avoid "bad" starting vector.
  arpackRestart<RC>(opt, "arpackmm.resid.out", nbDim, resid, false);
  a_int ldv = nbDim;
  RC * v = new RC[ldv*opt.nbCV]; for (a_int n = 0; n < ldv*opt.nbCV; n++) v[n] = zero; // Avoid "bad" starting vector.
  arpackRestart<RC>(opt, "arpackmm.v.out", ldv, v, true);
  a_int iparam[11];
  iparam[0] = 1; // Use exact shifts (=> we'll never have ido == 3).
  iparam[2] = opt.maxIt; // Maximum number of iterations.
  iparam[3] = 1; // Block size.
  iparam[4] = 0; // Number of ev found by arpack.
  iparam[6] = out.mode;
  a_int ipntr[14];
  RC * workd = new RC[3*nbDim]; for (a_int n = 0; n < 3*nbDim; n++) workd[n] = zero; // Avoid "bad" X/Y vector.
  a_int lworkl = opt.symPb ? opt.nbCV*opt.nbCV + 8*opt.nbCV : 3*opt.nbCV*opt.nbCV + 6*opt.nbCV;
  RC * workl = new RC[lworkl];
  a_int info = 0; // Use random initial residual vector.
  if (opt.restart) info = 1;

  // Initialize solver.

  auto start = chrono::high_resolution_clock::now();
  int rc = arpackMode<SLV, RC, FD>(opt, out.mode, A, B, solver);
  if (rc != 0) {cerr << "Error: bad arpack mode" << endl; return rc;}
  auto stop = chrono::high_resolution_clock::now();
  out.imsTime = chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;

  // Arpack solve.

  FD * rwork = NULL;
  do {
    // Call arpack.

    arpackAUPD(opt, &ido, bMat, nbDim, which, resid, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, &info);
    if (info ==  1) cerr << "Error: [dz][sn]aupd - KO: maximum number of iterations taken. Increase --maxIt..." << endl;
    if (info ==  3) cerr << "Error: [dz][sn]aupd - KO: no shifts could be applied. Increase --nbCV..." << endl;
    if (info == -9) cerr << "Error: [dz][sn]aupd - KO: starting vector is zero. Retry: play with shift..." << endl;
    if (info < 0) {cerr << "Error: [dz][sn]aupd - KO with info " << info << ", nbIt " << iparam[2] << endl; return 1;}

    // Reverse Communication Interface: perform actions according to arpack.

    start = chrono::high_resolution_clock::now();

    a_int xIdx = ipntr[0] - 1; // 0-based (Fortran is 1-based).
    a_int yIdx = ipntr[1] - 1; // 0-based (Fortran is 1-based).

    EV X(workd + xIdx, nbDim); // Arpack provides X.
    EV Y(workd + yIdx, nbDim); // Arpack provides Y.

    if (ido == -1) {
      if (iparam[6] == 1) {
        Y = A * X;
      }
      else if (iparam[6] == 2) {
        Y = A * X;
        auto YY = Y;          // Use copy of Y (not Y) for solve (avoid potential memory overwrite as Y is both in/out).
        Y = solver.solve(YY); // Y = B^-1 * A * X.
      }
      else if (iparam[6] == 3) {
        auto Z = B * X;      // Z = B * X.
        Y = solver.solve(Z); // Y = (A - sigma * B)^-1 * B * X.
      }
    }
    else if (ido == 1) {
      if (iparam[6] == 1) {
        Y = A * X;
      }
      else if (iparam[6] == 2) {
        Y = A * X;
        if (opt.symPb) X = Y; // Remark 5 in dsaupd documentation.
        auto YY = Y;          // Use copy of Y (not Y) for solve (avoid potential memory overwrite as Y is both in/out).
        Y = solver.solve(YY); // Y = B^-1 * A * X.
      }
      else if (iparam[6] == 3) {
        a_int zIdx = ipntr[2] - 1; // 0-based (Fortran is 1-based).
        EV Z(workd + zIdx, nbDim); // Arpack provides Z.
        Y = solver.solve(Z);       // Y = (A - sigma * B)^-1 * B * X.
      }
    }
    else if (ido == 2) {
      if      (iparam[6] == 1) Y =     X; // Y = I * X.
      else if (iparam[6] == 2) Y = B * X; // Y = B * X.
      else if (iparam[6] == 3) Y = B * X; // Y = B * X.
    }
    else if (ido != 99) {cerr << "Error: unexpected ido " << ido << " - KO" << endl; return 1;}

    stop = chrono::high_resolution_clock::now();
    out.rciTime += chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;

  } while (ido != 99);

  // Get arpack results (computed eigen values and vectors).

  out.nbIt = iparam[2]; // Actual number of iterations.
  bool rvec = true;
  char const * howmnyA = "A"; // Ritz vectors.
  char const * howmnyP = "P"; // Schur vectors.
  char const * howmny = opt.schur ? howmnyP : howmnyA;
  a_int * select = new a_int[opt.nbCV]; for (a_int n = 0; n < opt.nbCV; n++) select[n] = 1;
  a_int const nbZ = nbDim*(opt.nbEV+1); // Caution: opt.nbEV+1 for dneupd.
  RC * z = new RC[nbZ]; for (a_int n = 0; n < nbZ; n++) z[n] = zero;
  a_int ldz = nbDim;
  rc = arpackEUPD(opt, out, rvec, howmny, select, z, ldz, bMat, nbDim, which, resid, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info);
  if (rc != 0) {cerr << "Error: bad arpack eupd" << endl; return rc;}

  ofstream rfs("arpackmm.resid.out"); for (a_int n = 0; n < nbDim; n++) rfs << resid[n] << endl;
  ofstream vfs("arpackmm.v.out"); vfs << opt.nbCV << endl; for (a_int n = 0; n < ldv*opt.nbCV; n++) vfs << v[n] << endl;

  // Clean.

  if (rwork)  {delete [] rwork;  rwork  = NULL;}
  if (z)      {delete [] z;      z      = NULL;}
  if (select) {delete [] select; select = NULL;}
  if (workl)  {delete [] workl;  workl  = NULL;}
  if (workd)  {delete [] workd;  workd  = NULL;}
  if (v)      {delete [] v;      v      = NULL;}
  if (resid)  {delete [] resid;  resid  = NULL;}

  return 0;
}

template<typename EM>
int checkArpackEigVec(options const & opt, EM const & A, EM const & B, arpackEV const & out) {
  // Check eigen vectors.

  string rs = opt.schur ? "Schur" : "Ritz";

  if (opt.check && out.vec.size() == 0) {
    cerr << "Error: no " << rs << " value / vector found" << endl;
    return 1;
  }

  for (size_t i = 0; i < out.vec.size(); i++) {
    EigVecZ V = out.vec[i];
    complex<double> lambda = out.val[i];
    if (opt.verbose >= 1) {
      cout << endl;
      cout << rs << " value " << setw(3) << i << ": " << lambda << endl;
      if (opt.verbose >= 2) {
        cout << endl;
        cout << rs << " vector " << setw(3) << i << " (norm " << V.norm() << "): " << endl;
        cout << endl << V << endl;
      }
    }

    if (opt.check) {
      EigVecZ left = A.template cast<complex<double>>() * V;
      EigVecZ right = opt.stdPb ? V : B.template cast<complex<double>>() * V;
      right *= lambda;
      EigVecZ diff = left - right;
      if (diff.norm() > sqrt(opt.tol)) {
        cerr << endl << "Error: bad vector " << setw(3) << i << " (norm " << V.norm() << "):" << endl;
        cerr << endl << V << endl;
        cerr << endl << "Error: left side (A*V - norm " << left.norm() << "):" << endl;
        cerr << endl << left << endl;
        cerr << endl << "Error: right side (lambda*" << (opt.stdPb ? "" : "B*") << "V - norm " << right.norm() << "):" << endl;
        cerr << endl << right << endl;
        cerr << endl << "Error: diff (norm " << diff.norm() << ", sqrt(tol) " << sqrt(opt.tol) << "):" << endl;
        cerr << endl << diff << endl;
        return 1;
      }
      else {
        if (opt.verbose >= 1) {
          cout << endl << rs << " value/vector " << setw(3) << i << ": check OK";
          cout << ", diff (norm " << diff.norm() << ", sqrt(tol) " << sqrt(opt.tol) << ")" << endl;
        }
      }
    }
  }

  return 0;
}

void makeSigma(options const & opt,          float  & sigma) {sigma = (float) opt.sigmaReal;}

void makeSigma(options const & opt,         double  & sigma) {sigma = opt.sigmaReal;}

void makeSigma(options const & opt, complex< float> & sigma) {sigma = complex< float>((float) opt.sigmaReal, (float) opt.sigmaImag);}

void makeSigma(options const & opt, complex<double> & sigma) {sigma = complex<double>(opt.sigmaReal, opt.sigmaImag);}

template<typename RC, typename FD,
         typename EM, typename EV,
         typename SLVSLU, typename SLVSQR, typename SLVSLLT, typename SLVSLDLT>
int arpackSolveDrt(options & opt,
                   Eigen::SparseMatrix<RC> & A,
                   Eigen::SparseMatrix<RC> & B,
                   arpackEV & out) {
  // Solve with arpack using sparse matrices and direct solvers.

  int rc = 0;

  if (opt.slv.find("LU") != string::npos || opt.slv.find("QR") != string::npos) {
    stringstream clo(opt.slv);
    string slv; getline(clo, slv, '#');

    unique_ptr<double> slvDrtPvtThd = nullptr;
    string pivot; getline(clo, pivot);
    double pivotThd = 0.; stringstream pt(pivot); pt >> pivotThd;
    if (pt) { // Valid value read.
      slvDrtPvtThd = unique_ptr<double>(new double);
      if (slvDrtPvtThd) *slvDrtPvtThd = pivotThd;
    }

    if (slv == "LU") {
      SLVSLU solver;
      if (slvDrtPvtThd) solver.setPivotThreshold(*slvDrtPvtThd);
      rc = arpackSolve<RC, FD, EM, EV, SLVSLU>(opt, A, B, solver, out);
    }
    else if (slv == "QR") {
      SLVSQR solver;
      if (slvDrtPvtThd) solver.setPivotThreshold(*slvDrtPvtThd);
      rc = arpackSolve<RC, FD, EM, EV, SLVSQR>(opt, A, B, solver, out);
    }
    else {cerr << "Error: unknown solver - KO" << endl; return 1;}
  }
  else if (opt.slv.find("LLT") != string::npos || opt.slv.find("LDLT") != string::npos) {
    stringstream clo(opt.slv);
    string slv; getline(clo, slv, '#');

    unique_ptr<double> slvOffset = unique_ptr<double>(new double);
    string offset; getline(clo, offset, '#');
    double shiftOffset = 0.; stringstream so(offset); so >> shiftOffset; if (!so) shiftOffset = 0.;
    if (slvOffset) *slvOffset = shiftOffset;

    unique_ptr<double> slvScale = unique_ptr<double>(new double);
    string scale; getline(clo, scale);
    double shiftScale = 1.; stringstream ss(scale); ss >> shiftScale; if (!ss) shiftScale = 1.;
    if (slvScale) *slvScale = shiftScale;

    if (slv == "LLT") {
      SLVSLLT solver;
      if (slvOffset && slvScale) solver.setShift(*slvOffset, *slvScale);
      rc = arpackSolve<RC, FD, EM, EV, SLVSLLT >(opt, A, B, solver, out);
    }
    else if (slv == "LDLT") {
      SLVSLDLT solver;
      if (slvOffset && slvScale) solver.setShift(*slvOffset, *slvScale);
      rc = arpackSolve<RC, FD, EM, EV, SLVSLDLT>(opt, A, B, solver, out);
    }
    else {cerr << "Error: unknown solver - KO" << endl; return 1;}
  }
  else {cerr << "Error: unknown solver - KO" << endl; return 1;}

  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

  return 0;
}

template<typename RC, typename FD,
         typename EM, typename EV,
         typename SLVBCG, typename SLVBCGILU, typename SLVCG, typename SLVCGILU>
int arpackSolveItr(options & opt,
                   Eigen::SparseMatrix<RC> & A,
                   Eigen::SparseMatrix<RC> & B,
                   arpackEV & out) {
  // Solve with arpack using sparse matrices and iterative solvers.

  int rc = 0;

  if (opt.slv.find("BiCG") != string::npos || opt.slv.find("CG") != string::npos) {
    stringstream clo(opt.slvItrPC);
    string slvItrPC; getline(clo, slvItrPC, '#');

    unique_ptr<double> slvItrILUDropTol;
    if (slvItrPC == "ILU") {
      string dropTol; getline(clo, dropTol, '#');
      double iluDropTol = 1.; stringstream dt(dropTol); dt >> iluDropTol;
      if (dt) { // Valid value read.
        slvItrILUDropTol = unique_ptr<double>(new double);
        if (slvItrILUDropTol) *slvItrILUDropTol = iluDropTol;
      }
    }

    unique_ptr<int> slvItrILUFillFactor;
    if (slvItrPC == "ILU") {
      string fillFactor; getline(clo, fillFactor);
      int iluFillFactor = 2; stringstream ff(fillFactor); ff >> iluFillFactor;
      if (ff) { // Valid value read.
        slvItrILUFillFactor = unique_ptr<int>(new int);
        if (slvItrILUFillFactor) *slvItrILUFillFactor = iluFillFactor;
      }
    }

    if (opt.slv == "BiCG") {
      if (slvItrPC == "diag") {
        SLVBCG solver;
        if (opt.slvItrTol)   solver.setTolerance(*opt.slvItrTol);
        if (opt.slvItrMaxIt) solver.setMaxIterations(*opt.slvItrMaxIt);
        rc = arpackSolve<RC, FD, EM, EV, SLVBCG   >(opt, A, B, solver, out);
      }
      else if (slvItrPC == "ILU") {
        SLVBCGILU solver;
        if (opt.slvItrTol)       solver.setTolerance(*opt.slvItrTol);
        if (opt.slvItrMaxIt)     solver.setMaxIterations(*opt.slvItrMaxIt);
        if (slvItrILUDropTol)    solver.preconditioner().setDroptol(*slvItrILUDropTol);
        if (slvItrILUFillFactor) solver.preconditioner().setFillfactor(*slvItrILUFillFactor);
        rc = arpackSolve<RC, FD, EM, EV, SLVBCGILU>(opt, A, B, solver, out);
      }
      else {cerr << "Error: unknown preconditioner - KO" << endl; return 1;}
    }
    else if (opt.slv == "CG") {
      if (slvItrPC == "diag") {
        SLVCG solver;
        if (opt.slvItrTol)   solver.setTolerance(*opt.slvItrTol);
        if (opt.slvItrMaxIt) solver.setMaxIterations(*opt.slvItrMaxIt);
        rc = arpackSolve<RC, FD, EM, EV, SLVCG   >(opt, A, B, solver, out);
      }
      else if (slvItrPC == "ILU") {
        SLVCGILU solver;
        if (opt.slvItrTol)       solver.setTolerance(*opt.slvItrTol);
        if (opt.slvItrMaxIt)     solver.setMaxIterations(*opt.slvItrMaxIt);
        if (slvItrILUDropTol)    solver.preconditioner().setDroptol(*slvItrILUDropTol);
        if (slvItrILUFillFactor) solver.preconditioner().setFillfactor(*slvItrILUFillFactor);
        rc = arpackSolve<RC, FD, EM, EV, SLVCGILU>(opt, A, B, solver, out);
      }
      else {cerr << "Error: unknown preconditioner - KO" << endl; return 1;}
    }
  }
  else {cerr << "Error: unknown solver - KO" << endl; return 1;}

  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

  return 0;
}

template<typename RC, typename SLV>
void arpackSolveDrtSetThreshold(SLV & solver, unique_ptr<double> const & slvDrtPvtThd) {
  solver.setThreshold(Eigen::Default);
  if (slvDrtPvtThd) solver.setThreshold(*slvDrtPvtThd);
}

template<typename RC, typename SLV>
void arpackSolveDrtSetThreshold(Eigen::PartialPivLU<Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic>> & solver,
                                unique_ptr<double> const & slvDrtPvtThd) {
  // Thresholds only make sense for rank-revealing decompositions.
}

template<typename RC, typename SLV>
void arpackSolveDrtSetThreshold(Eigen::HouseholderQR<Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic>> & solver,
                                unique_ptr<double> const & slvDrtPvtThd) {
  // Thresholds only make sense for rank-revealing decompositions.
}

template<typename RC, typename FD,
         typename EM, typename EV,
         typename SLVDLU, typename SLVDQR, typename SLVDLLT, typename SLVDLDLT>
int arpackSolveDrt(options & opt,
                   Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> & A,
                   Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> & B,
                   arpackEV & out) {
  // Solve with arpack using dense matrices and direct solvers.

  int rc = 0;

  if (opt.slv.find("LU") != string::npos || opt.slv.find("QR") != string::npos) {
    stringstream clo(opt.slv);
    string slv; getline(clo, slv, '#');

    unique_ptr<double> slvDrtPvtThd = nullptr;
    string pivot; getline(clo, pivot);
    double pivotThd = 0.; stringstream pt(pivot); pt >> pivotThd;
    if (pt) { // Valid value read.
      slvDrtPvtThd = unique_ptr<double>(new double);
      if (slvDrtPvtThd) *slvDrtPvtThd = pivotThd;
    }

    a_int n = (out.mode == 1) ? 2 /*solver not used: reduce memory footprint*/ : A.rows();
    Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> S(n, n); S.setIdentity();
    if (slv == "LU") {
      SLVDLU solver(S);
      arpackSolveDrtSetThreshold<RC, SLVDLU>(solver, slvDrtPvtThd);
      rc = arpackSolve<RC, FD, EM, EV, SLVDLU>(opt, A, B, solver, out);
    }
    else if (slv == "QR") {
      SLVDQR solver(S);
      arpackSolveDrtSetThreshold<RC, SLVDQR>(solver, slvDrtPvtThd);
      rc = arpackSolve<RC, FD, EM, EV, SLVDQR>(opt, A, B, solver, out);
    }
    else {cerr << "Error: unknown solver - KO" << endl; return 1;}
  }
  else if (opt.slv.find("LLT") != string::npos || opt.slv.find("LDLT") != string::npos) {
    stringstream clo(opt.slv);
    string slv; getline(clo, slv, '#');

    a_int n = (out.mode == 1) ? 2 /*solver not used: reduce memory footprint*/ : A.rows();
    Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> S(n, n); S.setIdentity();
    if (slv == "LLT") {
      SLVDLLT solver(S);
      rc = arpackSolve<RC, FD, EM, EV, SLVDLLT >(opt, A, B, solver, out);
    }
    else if (slv == "LDLT") {
      SLVDLDLT solver(S);
      rc = arpackSolve<RC, FD, EM, EV, SLVDLDLT>(opt, A, B, solver, out);
    }
    else {cerr << "Error: unknown solver - KO" << endl; return 1;}
  }
  else {cerr << "Error: unknown solver - KO" << endl; return 1;}

  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

  return 0;
}

template<typename RC, typename FD,
         typename EM, typename ET, typename EV,
         typename SLV1, typename SLV2, typename SLV3, typename SLV4>
int arpackSolve(options & opt, arpackEV & out,
                int (*pf) (options &, EM & A, EM & B, arpackEV &) /* Pointer on (instance of) templated function */) {
  int rc = 0;

  // Read A and B matrices.

  EM A;
  auto start = chrono::high_resolution_clock::now();
  rc = createMatrix<RC, ET>(opt.fileA, A, opt.verbose, "A:");
  if (rc != 0) {cerr << "Error: read A KO" << endl; return rc;}
  auto stop = chrono::high_resolution_clock::now();
  double readATime = chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;
  cout << endl;
  cout << "INP: create A " << readATime << " s" << endl;

  if (opt.nbCV > A.cols()) opt.nbCV = A.cols(); /* Cut-off */

  EM B;
  if (!opt.stdPb) {
    start = chrono::high_resolution_clock::now();
    rc = createMatrix<RC, ET>(opt.fileB, B, opt.verbose, "B:");
    if (rc != 0) {cerr << "Error: read B KO" << endl; return rc;}
    stop = chrono::high_resolution_clock::now();
    double readBTime = chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;
    cout << endl;
    cout << "INP: create B " << readBTime << " s" << endl;

    if (A.rows() != B.rows()) {cerr << "Error: A.rows() != B.rows()" << endl; return rc;}
    if (A.cols() != B.cols()) {cerr << "Error: A.cols() != B.cols()" << endl; return rc;}
  }

  // If needed, transform the initial problem into a new one that arpack can handle.

  auto eps = numeric_limits<double>::epsilon();
  bool shiftReal = (opt.shiftReal && fabs(opt.sigmaReal) > eps) ? true : false;
  bool shiftImag = (opt.shiftImag && fabs(opt.sigmaImag) > eps) ? true : false;
  bool backTransform = false;
  int mode = 0;
  if (opt.stdPb) {
    mode = 1;
    if (shiftReal && !shiftImag) {
      EM I(A.rows(), A.cols());
      I.setIdentity();
      RC sigma; makeSigma(opt, sigma);
      A -= sigma*I;
      backTransform = true;
    }
  }
  else {
    mode = 2;
    if (shiftReal || shiftImag) mode = 3;
  }
  out.mode = mode;
  if (opt.verbose >= 1) {
    cout << endl;
    cout << "ARP: mode " << mode;
    cout << ", backTransform " << (backTransform ? "yes" : "no") << endl;
  }

  // Solve with arpack.
  // Note: when initializing Eigen solvers, API differ depending on solvers.

  rc = pf ? pf(opt, A, B, out) : 1; // Use function pointer to point at correct solver initialization.
  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

  // If needed, transform back the arpack problem into the initial problem.

  if (backTransform) {
    EM I(A.rows(), A.cols());
    I.setIdentity();
    RC sigma; makeSigma(opt, sigma);
    A += sigma*I;
    for (size_t i = 0; i < out.val.size(); i++) out.val[i] += sigma;
  }

  // Check.

  return checkArpackEigVec<EM>(opt, A, B, out);
}

int main(int argc, char ** argv) {
  // Check for options.

  options opt;
  int rc = opt.readCmdLine(argc, argv);
  if (rc != 0) {cerr << "Error: read cmd line KO" << endl; return rc;}
  cout << opt; // Print options.

  // Solve with arpack.

  sstats_c(); // Reset timers.
  sstatn_c(); // Reset timers.
  cstatn_c(); // Reset timers.

  bool itrSlv = true; // Use iterative solvers.
  if (opt.slv.find("LU")  != string::npos || opt.slv.find("QR")   != string::npos ||
      opt.slv.find("LLT") != string::npos || opt.slv.find("LDLT") != string::npos ) itrSlv = false;

  arpackEV out;
  auto start = chrono::high_resolution_clock::now();

  if (opt.dense) {
    if (itrSlv) {
      cerr << "Error: dense matrices does not support iterative solvers" << endl;
      return 1;
    }

    if (opt.simplePrec) {
      if (opt.cpxPb) {
        int (*pf) (options &, EigDMxC &, EigDMxC &, arpackEV &);
        if (opt.denseRR) {
          pf = arpackSolveDrt<complex< float>,  float, EigDMxC,          EigMpVC, EigDFLUC, EigDFQRC, EigDLLTC, EigDLDLTC>;
          rc = arpackSolve   <complex< float>,  float, EigDMxC, EigSMTC, EigMpVC, EigDFLUC, EigDFQRC, EigDLLTC, EigDLDLTC>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<complex< float>,  float, EigDMxC,          EigMpVC, EigDPLUC, EigDPQRC, EigDLLTC, EigDLDLTC>;
          rc = arpackSolve   <complex< float>,  float, EigDMxC, EigSMTC, EigMpVC, EigDPLUC, EigDPQRC, EigDLLTC, EigDLDLTC>(opt, out, pf);
        }
      }
      else {
        int (*pf) (options &, EigDMxS &, EigDMxS &, arpackEV &);
        if (opt.denseRR) {
          pf = arpackSolveDrt<         float ,  float, EigDMxS,          EigMpVS, EigDFLUS, EigDFQRS, EigDLLTS, EigDLDLTS>;
          rc = arpackSolve   <         float ,  float, EigDMxS, EigSMTS, EigMpVS, EigDFLUS, EigDFQRS, EigDLLTS, EigDLDLTS>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<         float ,  float, EigDMxS,          EigMpVS, EigDPLUS, EigDPQRS, EigDLLTS, EigDLDLTS>;
          rc = arpackSolve   <         float ,  float, EigDMxS, EigSMTS, EigMpVS, EigDPLUS, EigDPQRS, EigDLLTS, EigDLDLTS>(opt, out, pf);
        }
      }
    }
    else {
      if (opt.cpxPb) {
        int (*pf) (options &, EigDMxZ &, EigDMxZ &, arpackEV &);
        if (opt.denseRR) {
          pf = arpackSolveDrt<complex<double>, double, EigDMxZ,          EigMpVZ, EigDFLUZ, EigDFQRZ, EigDLLTZ, EigDLDLTZ>;
          rc = arpackSolve   <complex<double>, double, EigDMxZ, EigSMTZ, EigMpVZ, EigDFLUZ, EigDFQRZ, EigDLLTZ, EigDLDLTZ>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<complex<double>, double, EigDMxZ,          EigMpVZ, EigDPLUZ, EigDPQRZ, EigDLLTZ, EigDLDLTZ>;
          rc = arpackSolve   <complex<double>, double, EigDMxZ, EigSMTZ, EigMpVZ, EigDPLUZ, EigDPQRZ, EigDLLTZ, EigDLDLTZ>(opt, out, pf);
        }
      }
      else {
        int (*pf) (options &, EigDMxD &, EigDMxD &, arpackEV &);
        if (opt.denseRR) {
          pf = arpackSolveDrt<        double , double, EigDMxD,          EigMpVD, EigDFLUD, EigDFQRD, EigDLLTD, EigDLDLTD>;
          rc = arpackSolve   <        double , double, EigDMxD, EigSMTD, EigMpVD, EigDFLUD, EigDFQRD, EigDLLTD, EigDLDLTD>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<        double , double, EigDMxD,          EigMpVD, EigDPLUD, EigDPQRD, EigDLLTD, EigDLDLTD>;
          rc = arpackSolve   <        double , double, EigDMxD, EigSMTD, EigMpVD, EigDPLUD, EigDPQRD, EigDLLTD, EigDLDLTD>(opt, out, pf);
        }
      }
    }
  }
  else {
    if (opt.simplePrec) {
      if (opt.cpxPb) {
        int (*pf) (options &, EigSMxC &, EigSMxC &, arpackEV &);
        if (itrSlv) {
          pf = arpackSolveItr<complex< float>,  float, EigSMxC,          EigMpVC, EigSBiCGC, EigSBiCGILUC, EigSCGC,  EigSCGILUC>;
          rc = arpackSolve   <complex< float>,  float, EigSMxC, EigSMTC, EigMpVC, EigSBiCGC, EigSBiCGILUC, EigSCGC,  EigSCGILUC>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<complex< float>,  float, EigSMxC,          EigMpVC, EigSLUC,   EigSQRC,      EigSLLTC, EigSLDLTC >;
          rc = arpackSolve   <complex< float>,  float, EigSMxC, EigSMTC, EigMpVC, EigSLUC,   EigSQRC,      EigSLLTC, EigSLDLTC >(opt, out, pf);
        }
      }
      else {
        int (*pf) (options &, EigSMxS &, EigSMxS &, arpackEV &);
        if (itrSlv) {
          pf = arpackSolveItr<         float ,  float, EigSMxS,          EigMpVS, EigSBiCGS, EigSBiCGILUS, EigSCGS,  EigSCGILUS>;
          rc = arpackSolve   <         float ,  float, EigSMxS, EigSMTS, EigMpVS, EigSBiCGS, EigSBiCGILUS, EigSCGS,  EigSCGILUS>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<         float ,  float, EigSMxS,          EigMpVS, EigSLUS,   EigSQRS,      EigSLLTS, EigSLDLTS >;
          rc = arpackSolve   <         float ,  float, EigSMxS, EigSMTS, EigMpVS, EigSLUS,   EigSQRS,      EigSLLTS, EigSLDLTS >(opt, out, pf);
        }
      }
    }
    else {
      if (opt.cpxPb) {
        int (*pf) (options &, EigSMxZ &, EigSMxZ &, arpackEV &);
        if (itrSlv) {
          pf = arpackSolveItr<complex<double>, double, EigSMxZ,          EigMpVZ, EigSBiCGZ, EigSBiCGILUZ, EigSCGZ,  EigSCGILUZ>;
          rc = arpackSolve   <complex<double>, double, EigSMxZ, EigSMTZ, EigMpVZ, EigSBiCGZ, EigSBiCGILUZ, EigSCGZ,  EigSCGILUZ>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<complex<double>, double, EigSMxZ,          EigMpVZ, EigSLUZ,   EigSQRZ,      EigSLLTZ, EigSLDLTZ >;
          rc = arpackSolve   <complex<double>, double, EigSMxZ, EigSMTZ, EigMpVZ, EigSLUZ,   EigSQRZ,      EigSLLTZ, EigSLDLTZ >(opt, out, pf);
        }
      }
      else {
        int (*pf) (options &, EigSMxD &, EigSMxD &, arpackEV &);
        if (itrSlv) {
          pf = arpackSolveItr<        double , double, EigSMxD,          EigMpVD, EigSBiCGD, EigSBiCGILUD, EigSCGD,  EigSCGILUD>;
          rc = arpackSolve   <        double , double, EigSMxD, EigSMTD, EigMpVD, EigSBiCGD, EigSBiCGILUD, EigSCGD,  EigSCGILUD>(opt, out, pf);
        }
        else {
          pf = arpackSolveDrt<        double , double, EigSMxD,          EigMpVD, EigSLUD,   EigSQRD,      EigSLLTD, EigSLDLTD >;
          rc = arpackSolve   <        double , double, EigSMxD, EigSMTD, EigMpVD, EigSLUD,   EigSQRD,      EigSLLTD, EigSLDLTD >(opt, out, pf);
        }
      }
    }
  }
  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

  auto stop = chrono::high_resolution_clock::now();
  double fullTime = chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;
  cout << endl;
  cout << "OUT: mode " << out.mode << ", nb EV found " << out.val.size() << ", nb iterations " << out.nbIt << endl;
  cout << "OUT: init mode solver " << out.imsTime << " s, RCI time " << out.rciTime << " s" << endl;
  cout << "OUT: full time " << fullTime << " s" << endl;

  a_int nopx = 0, nbx = 0, nrorth = 0, nitref = 0, nrstrt = 0;
  float tsaupd = 0., tsaup2 = 0., tsaitr = 0., tseigt = 0., tsgets = 0., tsapps = 0., tsconv = 0.;
  float tnaupd = 0., tnaup2 = 0., tnaitr = 0., tneigt = 0., tngets = 0., tnapps = 0., tnconv = 0.;
  float tcaupd = 0., tcaup2 = 0., tcaitr = 0., tceigt = 0., tcgets = 0., tcapps = 0., tcconv = 0.;
  float tmvopx = 0., tmvbx = 0., tgetv0 = 0., titref = 0., trvec = 0.;
  stat_c(nopx, nbx, nrorth, nitref, nrstrt, tsaupd, tsaup2,
         tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd, tnaup2,
         tnaitr, tneigt, tngets, tnapps, tnconv, tcaupd, tcaup2,
         tcaitr, tceigt, tcgets, tcapps, tcconv, tmvopx, tmvbx,
         tgetv0, titref, trvec);
  cout << endl;
  cout << "STAT: total number of user OP*x operation                         " << nopx   << endl;
  cout << "STAT: total number of user  B*x operation                         " << nbx    << endl;
  cout << "STAT: total number of reorthogonalization steps taken             " << nrorth << endl;
  cout << "STAT: total number of it. refinement steps in reorthogonalization " << nitref << endl;
  cout << "STAT: total number of restart steps                               " << nrstrt << endl;

  return 0;
}

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=2 ts=2 et smartindent :*/
