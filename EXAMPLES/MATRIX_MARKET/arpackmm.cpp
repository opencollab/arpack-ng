#include <chrono>

#include "arpackSolver.hpp"
#include "debug_c.hpp"
#include "stat_c.hpp"

using namespace std;

class options {
 public:
  options() {
    fileA = "A.mtx";
    fileB = "N.A.";  // Not available.
    dense = false;
    denseRR = true;
    nbEV = 1;
    nbCV = 2 * nbEV + 1;
    stdPb = true;  // Standard or generalized (= not standard).
    symPb = true;
    cpxPb = false;
    simplePrec = false;  // Double precision.
    mag = string("LM");  // Large magnitude.
    shiftReal = false;
    shiftImag = false;
    sigmaReal = 0.;
    sigmaImag = 0.;  // Eigen value translation: look for lambda+sigma instead
                     // of lambda.
    invert =
        false;  // Eigen value invertion: look for 1./lambda instead of lambda.
    tol = 1.e-06;
    maxIt = 100;
    schur = false;  // Compute Ritz vectors.
    slv = "BiCG";
    slvItrTol = 1.e-6;
    slvItrMaxIt = 100;
    slvItrPC = "Diag";
    slvDrtPivot = 1.e-6;
    slvDrtOffset = 0.;
    slvDrtScale = 1.;
    check = true;
    verbose = 0;
    debug = 0;
    restart = false;
  };

  int readCmdLine(int argc, char** argv) {
    // Check for command line independent parameters.

    for (int a = 1; argv && a < argc; a++) {
      string clo = argv[a];  // Command line option.
      if (clo == "--help") return usage(0);
      if (clo == "--A") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        fileA = argv[a];
      }
      if (clo == "--dense") {
        dense = true;
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        string rr(argv[a]);
        if (rr != "true" && rr != "false") {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
        denseRR = (rr == "true") ? true : false;
      }
      if (clo == "--nbEV") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream nEV(argv[a]);
        nEV >> nbEV;
        if (!nEV) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
        nbCV = 2 * nbEV + 1;
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
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        mag = argv[a];  // small mag (likely poor perf) <=> large mag + invert
                        // (likely good perf).
        bool ok = (mag == "LM" || mag == "SM" || mag == "LR" || mag == "SR" ||
                   mag == "LI" || mag == "SI")
                      ? true
                      : false;
        if (!ok) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--shiftReal") {
        shiftReal = true;
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream s(argv[a]);
        s >> sigmaReal;
        if (!s) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--shiftImag") {
        shiftImag = true;
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream s(argv[a]);
        s >> sigmaImag;
        if (!s) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--invert") invert = true;
      if (clo == "--tol") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream t(argv[a]);
        t >> tol;
        if (!t) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--maxIt") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream mi(argv[a]);
        mi >> maxIt;
        if (!mi) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--slv") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        slv = argv[a];
      }
      if (clo == "--slvItrTol") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream t(argv[a]);
        double tol = 0.;
        t >> slvItrTol;
        if (!t) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--slvItrMaxIt") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream mi(argv[a]);
        int maxIt = 0;
        mi >> slvItrMaxIt;
        if (!mi) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--slvItrPC") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream pc(argv[a]);
        pc >> slvItrPC;
        if (!pc) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--slvDrtPivot") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream pv(argv[a]);
        pv >> slvDrtPivot;
        if (!pv) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--slvDrtOffset") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream of(argv[a]);
        of >> slvDrtOffset;
        if (!of) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--slvDrtScale") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream sc(argv[a]);
        sc >> slvDrtScale;
        if (!sc) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--noCheck") check = false;
      if (clo == "--verbose") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream vb(argv[a]);
        vb >> verbose;
        if (!vb) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--debug") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream dbg(argv[a]);
        dbg >> debug;
        if (!dbg) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
        if (debug > 3) debug = 3;
        debug_c(6, -6, debug, debug, debug, debug, debug, debug, debug, debug,
                debug, debug, debug, debug, debug, debug, debug, debug, debug,
                debug, debug, debug, debug, debug);
      }
      if (clo == "--restart") restart = true;
    }

    // Check for command line dependent parameters.

    for (int a = 1; argv && a < argc; a++) {
      string clo = argv[a];  // Command line option.
      if (clo == "--nbCV") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        stringstream nCV(argv[a]);
        nCV >> nbCV;
        if (!nCV) {
          cerr << "Error: bad " << clo << " - bad argument" << endl;
          return usage();
        }
      }
      if (clo == "--B") {
        a++;
        if (a >= argc) {
          cerr << "Error: bad " << clo << " - need argument" << endl;
          return usage();
        }
        fileB = argv[a];
      }
    }

    // Sanity checks.

    if (!stdPb && fileB.empty()) {
      cerr << "Error: generalized problem without B matrix" << endl;
      return usage();
    }

    return 0;
  };

  int usage(int rc = 1) {
    cout << "Usage: running arpack with matrix market files to check for eigen "
            "values/vectors."
         << endl;
    cout << endl;
    cout << "  --A F:            file name of matrix A such that A X = lambda "
            "X. (standard)"
         << endl;
    cout << "                    the file F must be compliant with the matrix "
            "market format."
         << endl;
    cout << "                    default: A.mtx" << endl;
    cout << "  --B F:            file name of matrix B such that A X = lambda "
            "B X. (generalized)"
         << endl;
    cout << "                    the file F must be compliant with the matrix "
            "market format."
         << endl;
    cout << "                    default: N.A. for standard problem, or, B.mtx "
            "for generalized problem"
         << endl;
    cout << "  --dense RR:       consider A and B as dense matrices." << endl;
    cout << "                    if RR = true,  use more-stable-but-slow "
            "versions of LU / QR (rank revealing)."
         << endl;
    cout << "                    if RR = false, use less-stable-but-fast "
            "versions of LU / QR (depends on condition number)."
         << endl;
    cout << "                    Notes:" << endl;
    cout << "                      - only direct solvers are available when "
            "using dense matrices."
         << endl;
    cout
        << "                      - RR does not impact the use of LLT and LDLT."
        << endl;
    cout << "                      - thresholds only make sense for "
            "rank-revealing decompositions."
         << endl;
    cout << "                    default: consider A and B as sparse matrices"
         << endl;
    cout << "  --nbEV:           number of eigen values/vectors to compute."
         << endl;
    cout << "                    default: 1" << endl;
    cout << "  --nbCV:           number of columns of the matrix V." << endl;
    cout << "                    default: 2*nbEV+1" << endl;
    cout << "  --genPb:          generalized problem." << endl;
    cout << "                    default: standard problem" << endl;
    cout << "  --nonSymPb:       non symmetric problem (<=> use dn[ae]upd)."
         << endl;
    cout << "                    default: symmetric problem (<=> use ds[ae]upd)"
         << endl;
    cout << "  --cpxPb:          complex (non symmetric) problem (<=> use "
            "zn[ae]upd)."
         << endl;
    cout << "                    default: false (<=> use d*[ae]upd)" << endl;
    cout << "  --simplePrec:     use simple precision (less accurate, but, "
            "half memory footprint)."
         << endl;
    cout << "                    default: false (<=> use double precision: use "
            "[dz]*upd)"
         << endl;
    cout << "  --mag M:          set magnitude of eigen values to look for "
            "(LM, SM, LR, SR, LI, SI)."
         << endl;
    cout << "                    default: large magnitude (LM)" << endl;
    cout << "  --shiftReal S:    real shift where sigma = S (look for lambda+S "
            "instead of lambda)."
         << endl;
    cout << "                    default: no shift, S = 0." << endl;
    cout << "  --shiftImag S:    imaginary shift where sigma = S (look for "
            "lambda+S instead of lambda)."
         << endl;
    cout << "                    default: no shift, S = 0." << endl;
    cout << "  --invert:         invert mode (look for 1./lambda instead of "
            "lambda)."
         << endl;
    cout << "                    default: no invert" << endl;
    cout << "  --tol T:          tolerance T." << endl;
    cout << "                    default: 1.e-06" << endl;
    cout << "  --maxIt M:        maximum iterations M." << endl;
    cout << "                    default: 100" << endl;
    cout << "  --schur:          compute Schur vectors." << endl;
    cout << "                      the Schur decomposition is such that A = "
            "Q^H x T x Q where:"
         << endl;
    cout << "                        - the H superscript refers to the "
            "Hermitian transpose: Q^H = (Q^t)^*."
         << endl;
    cout
        << "                        - Q is unitary: Q is such that Q^H x Q = I."
        << endl;
    cout << "                        - T is an upper-triangular matrix whose "
            "diagonal elements are the eigenvalues of A."
         << endl;
    cout << "                      every square matrix has a Schur "
            "decomposition: columns of Q are the Schur vectors."
         << endl;
    cout << "                      for a general matrix A, there is no "
            "relation between Schur vectors of A and eigenvectors of A."
         << endl;
    cout << "                      if q_j is the j-th Schur vector, then A x "
            "q_j is a linear combination of q_1, ..., q_j."
         << endl;
    cout << "                      Schur vectors q_1, q_2, ..., q_j span an "
            "invariant subspace of A."
         << endl;
    cout << "                      the Schur vectors and eigenvectors of A are "
            "the same if A is a normal matrix."
         << endl;
    cout << "                    default: compute Ritz vectors (approximations "
            "of eigen vectors)"
         << endl;
    cout << "  --slv S:          solver (needed if arpack mode > 1)." << endl;
    cout << "                      BiCG: iterative method, any matrices"
         << endl;
    cout << "                      CG:   iterative method, sym matrices only"
         << endl;
    cout << "                      LU:   direct method, any matrices (pivoting "
            "needed)"
         << endl;
    cout << "                      QR:   direct method, any matrices (pivoting "
            "needed)"
         << endl;
    cout << "                      LLT:  direct method, SPD matrices only "
            "(pivoting not needed)"
         << endl;
    cout << "                      LDLT: direct method, symmetric positive "
            "semi-definite matrices only (pivoting not needed)"
         << endl;
    cout << "                    default: BiCG" << endl;
    cout << "  --slvItrTol T:    solver tolerance T (for iterative solvers)."
         << endl;
    cout << "                    default: 1.e-6" << endl;
    cout << "  --slvItrMaxIt M:  solver maximum iterations M (for iterative "
            "solvers)."
         << endl;
    cout << "                    default: 100" << endl;
    cout << "  --slvItrPC PC:    solver preconditioner (for iterative solvers)."
         << endl;
    cout << "                      PC preconditioner:" << endl;
    cout << "                        Diag:    eigen diagonal preconditioner "
            "(Jacobi)."
         << endl;
    cout << "                        ILU#D#F: eigen ILU preconditioner."
         << endl;
    cout << "                          D:       drop tolerance." << endl;
    cout << "                          F:       fill factor." << endl;
    cout << "                    default: diagonal preconditioner (Jacobi)"
         << endl;
    cout << "  --slvDrtPivot: P  solver pivot P (for direct solvers)." << endl;
    cout << "                    default: 1.e-06" << endl;
    cout << "  --slvDrtOffset: O solver offset O (for direct solvers)." << endl;
    cout << "                    default: 0." << endl;
    cout << "  --slvDrtScale: S  solver scale S (for direct solvers)." << endl;
    cout << "                    default: 1." << endl;
    cout << "  --noCheck:        check arpack eigen values/vectors." << endl;
    cout << "                    check will fail if Schur vectors are computed "
            "and A is NOT a normal matrix."
         << endl;
    cout << "                    default: check" << endl;
    cout << "  --verbose V:      verbosity level (up to 3)." << endl;
    cout << "                    default: 0" << endl;
    cout << "  --debug D:        debug level (up to 3)." << endl;
    cout << "                    default: 0" << endl;
    cout << "  --restart:        restart from previous run (which had produced "
            "arpackSolver.*.out)."
         << endl;
    cout << "                    restart from eigen basis approximation "
            "computed during a previous run."
         << endl;
    cout << "                    default: false" << endl;
    if (rc == 0) exit(0);
    return rc;
  };

  friend ostream& operator<<(ostream& ostr, options const& opt);

  string fileA;
  string fileB;
  bool dense;
  bool denseRR;
  a_int nbEV;
  a_int nbCV;
  bool stdPb;  // Standard or generalized (= not standard).
  bool symPb;
  bool cpxPb;
  bool simplePrec;
  string mag;  // Magnitude <=> "which" arpack parameter.
  bool shiftReal, shiftImag;
  double sigmaReal, sigmaImag;  // Eigen value translation: look for
                                // lambda+sigma instead of lambda.
  bool invert;  // Eigen value invertion: look for 1./lambda instead of lambda.
  double tol;
  int maxIt;
  bool schur;
  string slv;
  double slvItrTol;
  int slvItrMaxIt;
  string slvItrPC;
  double slvDrtPivot;
  double slvDrtOffset;
  double slvDrtScale;
  bool check;
  int verbose;
  a_int debug;
  bool restart;
};

ostream& operator<<(ostream& ostr, options const& opt) {
  ostr << "OPT: A " << opt.fileA << ", B " << opt.fileB;
  if (opt.dense && opt.denseRR)
    ostr << ", dense yes (RR true)";
  else if (opt.dense && !opt.denseRR)
    ostr << ", dense yes (RR false)";
  else
    ostr << ", dense no";
  ostr << ", nbEV " << opt.nbEV << ", nbCV " << opt.nbCV << ", stdPb "
       << (opt.stdPb ? "yes" : "no");
  ostr << ", symPb " << (opt.symPb ? "yes" : "no") << ", cpxPb "
       << (opt.cpxPb ? "yes" : "no");
  ostr << ", simplePrec " << (opt.simplePrec ? "yes" : "no") << ", mag "
       << opt.mag << endl;
  ostr << "OPT: shiftReal " << (opt.shiftReal ? "yes" : "no") << ", sigmaReal "
       << opt.sigmaReal;
  ostr << ", shiftImag " << (opt.shiftImag ? "yes" : "no") << ", sigmaImag "
       << opt.sigmaImag;
  ostr << ", invert " << (opt.invert ? "yes" : "no") << ", tol " << opt.tol
       << ", maxIt " << opt.maxIt;
  ostr << ", " << (opt.schur ? "Schur" : "Ritz") << " vectors" << endl;
  ostr << "OPT: slv " << opt.slv << ", slvItrPC " << opt.slvItrPC
       << ", slvItrTol " << opt.slvItrTol;
  ostr << ", slvItrMaxIt " << opt.slvItrMaxIt << ", slvDrtPivot "
       << opt.slvDrtPivot;
  ostr << ", slvDrtOffset " << opt.slvDrtOffset << ", slvDrtScale "
       << opt.slvDrtScale << endl;
  ostr << "OPT: check " << (opt.check ? "yes" : "no") << ", verbose "
       << opt.verbose << ", debug " << opt.debug;
  ostr << ", restart " << (opt.restart ? "yes" : "no") << endl;
  return ostr;
}

class output {
 public:
  output() {
    nbVal = 0;
    mode = 0;
    nbIt = 0;
    imsTime = 0.;
    rciTime = 0.;
  };

  int nbVal;       // Eigen values.
  int mode;        // Arpack mode.
  int nbIt;        // Arpack number of iterations.
  double imsTime;  // Init mode solver time.
  double rciTime;  // Reverse communication interface time.
};

template <typename RC, typename FD, typename EM, typename SLV>
int itrSolve(options& opt, output& out, double const& slvItrILUDropTol,
             double const& slvItrILUFillFactor) {
  // Init solver.

  arpackItrSolver<RC, FD, EM, SLV> as;
  as.symPb = opt.symPb;
  as.nbEV = opt.nbEV;
  as.nbCV = opt.nbCV;
  as.tol = opt.tol;
  as.sigmaReal = opt.sigmaReal;
  as.sigmaImag = opt.sigmaImag;
  as.dumpToFile = true;
  as.restartFromFile = opt.restart;
  as.mag = opt.mag;
  as.maxIt = opt.maxIt;
  as.schur = opt.schur;
  as.verbose = opt.verbose;
  as.slvTol = opt.slvItrTol;
  as.slvMaxIt = opt.slvItrMaxIt;
  as.slvILUDropTol = slvItrILUDropTol;
  as.slvILUFillFactor = slvItrILUFillFactor;

  // Read A and B matrices.

  EM A;
  auto start = chrono::high_resolution_clock::now();
  int rc = as.createMatrix(opt.fileA, A);
  if (rc != 0) {
    cerr << "Error: read A KO" << endl;
    return rc;
  }
  auto stop = chrono::high_resolution_clock::now();
  double readATime =
      chrono::duration_cast<chrono::milliseconds>(stop - start).count() / 1000.;
  cout << endl;
  cout << "INP: create A " << readATime << " s" << endl;

  if (opt.nbCV > A.cols()) opt.nbCV = A.cols(); /* Cut-off */

  EM B;
  if (!opt.stdPb) {
    start = chrono::high_resolution_clock::now();
    rc = as.createMatrix(opt.fileB, B);
    if (rc != 0) {
      cerr << "Error: read B KO" << endl;
      return rc;
    }
    stop = chrono::high_resolution_clock::now();
    double readBTime =
        chrono::duration_cast<chrono::milliseconds>(stop - start).count() /
        1000.;
    cout << endl;
    cout << "INP: create B " << readBTime << " s" << endl;

    if (A.rows() != B.rows()) {
      cerr << "Error: A.rows() != B.rows()" << endl;
      return rc;
    }
    if (A.cols() != B.cols()) {
      cerr << "Error: A.cols() != B.cols()" << endl;
      return rc;
    }
  }

  // Solve.

  rc = as.solve(A, opt.stdPb ? nullptr : &B);
  if (rc != 0) {
    cerr << "Error: solve KO" << endl;
    return rc;
  }
  if (opt.check) {
    rc = as.checkEigVec(A, opt.stdPb ? nullptr : &B);
    if (rc != 0) {
      cerr << "Error: check KO" << endl;
      return rc;
    }
  }

  // Retrieve outputs.

  out.nbVal = as.val.size();
  out.mode = as.mode;
  out.nbIt = as.nbIt;
  out.imsTime = as.imsTime;
  out.rciTime = as.rciTime;

  return 0;
}

template <typename RC, typename FD, typename EM, typename SLV>
int drtSolve(options& opt, output& out) {
  // Init solver.

  arpackDrtSolver<RC, FD, EM, SLV> as;
  as.symPb = opt.symPb;
  as.nbEV = opt.nbEV;
  as.nbCV = opt.nbCV;
  as.tol = opt.tol;
  as.sigmaReal = opt.sigmaReal;
  as.sigmaImag = opt.sigmaImag;
  as.dumpToFile = true;
  as.restartFromFile = opt.restart;
  as.mag = opt.mag;
  as.maxIt = opt.maxIt;
  as.schur = opt.schur;
  as.verbose = opt.verbose;
  as.slvPvtThd = opt.slvDrtPivot;
  as.slvOffset = opt.slvDrtOffset;
  as.slvScale = opt.slvDrtScale;

  // Read A and B matrices.

  EM A;
  auto start = chrono::high_resolution_clock::now();
  int rc = as.createMatrix(opt.fileA, A);
  if (rc != 0) {
    cerr << "Error: read A KO" << endl;
    return rc;
  }
  auto stop = chrono::high_resolution_clock::now();
  double readATime =
      chrono::duration_cast<chrono::milliseconds>(stop - start).count() / 1000.;
  cout << endl;
  cout << "INP: create A " << readATime << " s" << endl;

  if (opt.nbCV > A.cols()) opt.nbCV = A.cols(); /* Cut-off */

  EM B;
  if (!opt.stdPb) {
    start = chrono::high_resolution_clock::now();
    rc = as.createMatrix(opt.fileB, B);
    if (rc != 0) {
      cerr << "Error: read B KO" << endl;
      return rc;
    }
    stop = chrono::high_resolution_clock::now();
    double readBTime =
        chrono::duration_cast<chrono::milliseconds>(stop - start).count() /
        1000.;
    cout << endl;
    cout << "INP: create B " << readBTime << " s" << endl;

    if (A.rows() != B.rows()) {
      cerr << "Error: A.rows() != B.rows()" << endl;
      return rc;
    }
    if (A.cols() != B.cols()) {
      cerr << "Error: A.cols() != B.cols()" << endl;
      return rc;
    }
  }

  // Solve.

  rc = as.solve(A, opt.stdPb ? nullptr : &B);
  if (rc != 0) {
    cerr << "Error: solve KO" << endl;
    return rc;
  }
  if (opt.check) {
    rc = as.checkEigVec(A, opt.stdPb ? nullptr : &B);
    if (rc != 0) {
      cerr << "Error: check KO" << endl;
      return rc;
    }
  }

  // Retrieve outputs.

  out.nbVal = as.val.size();
  out.mode = as.mode;
  out.nbIt = as.nbIt;
  out.imsTime = as.imsTime;
  out.rciTime = as.rciTime;

  return 0;
}

template <typename RC, typename FD, typename EM, typename SLV1, typename SLV2,
          typename SLV3, typename SLV4>
int drtSolve(options& opt, output& out) {
  int rc = 1;

  if (opt.slv == "LU") rc = drtSolve<RC, FD, EM, SLV1>(opt, out);
  if (opt.slv == "QR") rc = drtSolve<RC, FD, EM, SLV2>(opt, out);
  if (opt.slv == "LLT") rc = drtSolve<RC, FD, EM, SLV3>(opt, out);
  if (opt.slv == "LDLT") rc = drtSolve<RC, FD, EM, SLV4>(opt, out);

  return rc;
}

template <typename RC, typename FD, typename EM, typename SLV1, typename SLV2,
          typename SLV3, typename SLV4>
int itrSolve(options& opt, output& out) {
  int rc = 1;

  stringstream clo(opt.slvItrPC);
  string slvItrPC;
  getline(clo, slvItrPC, '#');

  double slvItrILUDropTol = 1.;
  if (slvItrPC == "ILU") {
    string dropTol;
    getline(clo, dropTol, '#');
    stringstream dt(dropTol);
    dt >> slvItrILUDropTol;
  }

  int slvItrILUFillFactor = 2;
  if (slvItrPC == "ILU") {
    string fillFactor;
    getline(clo, fillFactor);
    stringstream ff(fillFactor);
    ff >> slvItrILUFillFactor;
  }

  if (opt.slv == "BiCG") {
    if (slvItrPC == "Diag")
      rc = itrSolve<RC, FD, EM, SLV1>(opt, out, slvItrILUDropTol,
                                      slvItrILUFillFactor);
    if (slvItrPC == "ILU")
      rc = itrSolve<RC, FD, EM, SLV2>(opt, out, slvItrILUDropTol,
                                      slvItrILUFillFactor);
  }
  if (opt.slv == "CG") {
    if (slvItrPC == "Diag")
      rc = itrSolve<RC, FD, EM, SLV3>(opt, out, slvItrILUDropTol,
                                      slvItrILUFillFactor);
    if (slvItrPC == "ILU")
      rc = itrSolve<RC, FD, EM, SLV4>(opt, out, slvItrILUDropTol,
                                      slvItrILUFillFactor);
  }

  return rc;
}

int main(int argc, char** argv) {
  // Check for options.

  options opt;
  int rc = opt.readCmdLine(argc, argv);
  if (rc != 0) {
    cerr << "Error: read cmd line KO" << endl;
    return rc;
  }
  cout << opt;  // Print options.

  // Solve with arpack.

  sstats_c();  // Reset timers.
  sstatn_c();  // Reset timers.
  cstatn_c();  // Reset timers.

  bool itrSlv = true;  // Use iterative solvers.
  if (opt.slv.find("LU") != string::npos ||
      opt.slv.find("QR") != string::npos ||
      opt.slv.find("LLT") != string::npos ||
      opt.slv.find("LDLT") != string::npos)
    itrSlv = false;

  output out;
  auto start = chrono::high_resolution_clock::now();
  if (opt.dense) {
    if (itrSlv) {
      cerr << "Error: dense matrices does not support iterative solvers"
           << endl;
      return 1;
    }

    if (opt.simplePrec) {
      if (opt.cpxPb) {
        if (opt.denseRR) {
          rc = drtSolve<complex<float>, float, EigDMxC, EigDFLUC, EigDFQRC,
                        EigDLLTC, EigDLDLTC>(opt, out);
        } else {
          rc = drtSolve<complex<float>, float, EigDMxC, EigDPLUC, EigDPQRC,
                        EigDLLTC, EigDLDLTC>(opt, out);
        }
      } else {
        if (opt.denseRR) {
          rc = drtSolve<float, float, EigDMxS, EigDFLUS, EigDFQRS, EigDLLTS,
                        EigDLDLTS>(opt, out);
        } else {
          rc = drtSolve<float, float, EigDMxS, EigDPLUS, EigDPQRS, EigDLLTS,
                        EigDLDLTS>(opt, out);
        }
      }
    } else {
      if (opt.cpxPb) {
        if (opt.denseRR) {
          rc = drtSolve<complex<double>, double, EigDMxZ, EigDFLUZ, EigDFQRZ,
                        EigDLLTZ, EigDLDLTZ>(opt, out);
        } else {
          rc = drtSolve<complex<double>, double, EigDMxZ, EigDPLUZ, EigDPQRZ,
                        EigDLLTZ, EigDLDLTZ>(opt, out);
        }
      } else {
        if (opt.denseRR) {
          rc = drtSolve<double, double, EigDMxD, EigDFLUD, EigDFQRD, EigDLLTD,
                        EigDLDLTD>(opt, out);
        } else {
          rc = drtSolve<double, double, EigDMxD, EigDPLUD, EigDPQRD, EigDLLTD,
                        EigDLDLTD>(opt, out);
        }
      }
    }
  } else {
    if (opt.simplePrec) {
      if (opt.cpxPb) {
        if (itrSlv) {
          rc = itrSolve<complex<float>, float, EigSMxC, EigSBiCGC, EigSBiCGILUC,
                        EigSCGC, EigSCGILUC>(opt, out);
        } else {
          rc = drtSolve<complex<float>, float, EigSMxC, EigSLUC, EigSQRC,
                        EigSLLTC, EigSLDLTC>(opt, out);
        }
      } else {
        if (itrSlv) {
          rc = itrSolve<float, float, EigSMxS, EigSBiCGS, EigSBiCGILUS, EigSCGS,
                        EigSCGILUS>(opt, out);
        } else {
          rc = drtSolve<float, float, EigSMxS, EigSLUS, EigSQRS, EigSLLTS,
                        EigSLDLTS>(opt, out);
        }
      }
    } else {
      if (opt.cpxPb) {
        if (itrSlv) {
          rc = itrSolve<complex<double>, double, EigSMxZ, EigSBiCGZ,
                        EigSBiCGILUZ, EigSCGZ, EigSCGILUZ>(opt, out);
        } else {
          rc = drtSolve<complex<double>, double, EigSMxZ, EigSLUZ, EigSQRZ,
                        EigSLLTZ, EigSLDLTZ>(opt, out);
        }
      } else {
        if (itrSlv) {
          rc = itrSolve<double, double, EigSMxD, EigSBiCGD, EigSBiCGILUD,
                        EigSCGD, EigSCGILUD>(opt, out);
        } else {
          rc = drtSolve<double, double, EigSMxD, EigSLUD, EigSQRD, EigSLLTD,
                        EigSLDLTD>(opt, out);
        }
      }
    }
  }
  if (rc != 0) {
    cerr << "Error: arpack solve KO" << endl;
    return rc;
  }

  // Output results and stats.

  auto stop = chrono::high_resolution_clock::now();
  double fullTime =
      chrono::duration_cast<chrono::milliseconds>(stop - start).count() / 1000.;
  cout << endl;
  cout << "OUT: mode " << out.mode << ", nb EV found " << out.nbVal
       << ", nb iterations " << out.nbIt << endl;
  cout << "OUT: init mode solver " << out.imsTime << " s, RCI time "
       << out.rciTime << " s" << endl;
  cout << "OUT: full time " << fullTime << " s" << endl;

  a_int nopx = 0, nbx = 0, nrorth = 0, nitref = 0, nrstrt = 0;
  float tsaupd = 0., tsaup2 = 0., tsaitr = 0., tseigt = 0., tsgets = 0.,
        tsapps = 0., tsconv = 0.;
  float tnaupd = 0., tnaup2 = 0., tnaitr = 0., tneigt = 0., tngets = 0.,
        tnapps = 0., tnconv = 0.;
  float tcaupd = 0., tcaup2 = 0., tcaitr = 0., tceigt = 0., tcgets = 0.,
        tcapps = 0., tcconv = 0.;
  float tmvopx = 0., tmvbx = 0., tgetv0 = 0., titref = 0., trvec = 0.;
  stat_c(nopx, nbx, nrorth, nitref, nrstrt, tsaupd, tsaup2, tsaitr, tseigt,
         tsgets, tsapps, tsconv, tnaupd, tnaup2, tnaitr, tneigt, tngets, tnapps,
         tnconv, tcaupd, tcaup2, tcaitr, tceigt, tcgets, tcapps, tcconv, tmvopx,
         tmvbx, tgetv0, titref, trvec);
  cout << endl;
  cout << "STAT: total number of user OP*x operation                         "
       << nopx << endl;
  cout << "STAT: total number of user  B*x operation                         "
       << nbx << endl;
  cout << "STAT: total number of reorthogonalization steps taken             "
       << nrorth << endl;
  cout << "STAT: total number of it. refinement steps in reorthogonalization "
       << nitref << endl;
  cout << "STAT: total number of restart steps                               "
       << nrstrt << endl;

  return 0;
}

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=2 ts=2 et smartindent :*/
