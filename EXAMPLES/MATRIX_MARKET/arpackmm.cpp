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
#include "arpack.h"
#include "debug_c.hpp"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

using namespace std;

typedef Eigen::SparseMatrix<        double>                     EigMatR; // Real.
typedef Eigen::Triplet     <        double>                     EigCooR; // Real.
typedef Eigen::SparseMatrix<complex<double>>                    EigMatC; // Complex.
typedef Eigen::Triplet     <complex<double>>                    EigCooC; // Complex.
typedef Eigen::Matrix      <        double,  Eigen::Dynamic, 1> EigVecR; // Real.
typedef Eigen::Map         <EigVecR>                            EigMpVR; // Real.
typedef Eigen::Matrix      <complex<double>, Eigen::Dynamic, 1> EigVecC; // Complex.
typedef Eigen::Map         <EigVecC>                            EigMpVC; // Complex.
typedef Eigen::BiCGSTAB         <EigMatR>                       EigBiCG;
typedef Eigen::ConjugateGradient<EigMatR>                       EigCG;
typedef Eigen::SparseLU<EigMatR, Eigen::COLAMDOrdering<int>>    EigSLU;
typedef Eigen::SparseQR<EigMatR, Eigen::COLAMDOrdering<int>>    EigSQR;

class options {
  public:
    options() {
      fileA = "A.mtx";
      fileB = "N.A."; // Not available.
      nbEV = 1;
      nbCV = 2*nbEV + 1;
      stdPb = true; // Standard or generalized (= not standard).
      symPb = true;
      mag = string("LM"); // Large magnitude.
      shiftReal = false; shiftImag = false;
      sigmaReal = 0.; sigmaImag = 0.; // Eigen value translation: look for lambda+sigma instead of lambda.
      invert = false; // Eigen value invertion: look for 1./lambda instead of lambda.
      tol = 1.e-06;
      maxIt = 100;
      slv = "BiCG";
      slvTol = 1.e-06;
      slvMaxIt = 100;
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
        if (clo == "--slvTol") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream t(argv[a]);
          t >> slvTol; if (!t) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
        }
        if (clo == "--slvMaxIt") {
          a++; if (a >= argc) {cerr << "Error: bad " << clo << " - need argument" << endl; return usage();}
          stringstream mi(argv[a]);
          mi >> slvMaxIt; if (!mi) {cerr << "Error: bad " << clo << " - bad argument" << endl; return usage();}
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
      cout << "Usage: running arpack to check for eigen values/vectors." << endl;
      cout << endl;
      cout << "  --A F:         file name of matrix A such that A X = lambda X. (standard)" << endl;
      cout << "                 default: A.mtx" << endl;
      cout << "  --B F:         file name of matrix B such that A X = lambda B X. (generalized)" << endl;
      cout << "                 default: N.A. for standard problem, or, B.mtx for generalized problem" << endl;
      cout << "  --nbEV:        number of eigen values/vectors to compute." << endl;
      cout << "                 default: 1" << endl;
      cout << "  --nbCV:        number of columns of the matrix V." << endl;
      cout << "                 default: 2*nbEV+1" << endl;
      cout << "  --genPb:       generalized problem." << endl;
      cout << "                 default: standard problem" << endl;
      cout << "  --nonSymPb:    non symmetric problem." << endl;
      cout << "                 default: symmetric problem" << endl;
      cout << "  --mag M:       set magnitude of eigen values to look for (LM, SM, LR, SR, LI, SI)." << endl;
      cout << "                 default: large magnitude (LM)" << endl;
      cout << "  --shiftReal S: real shift where sigma = S (look for lambda+S instead of lambda)." << endl;
      cout << "                 default: no shift, S = 0." << endl;
      cout << "  --shiftImag S: imaginary shift where sigma = S (look for lambda+S instead of lambda)." << endl;
      cout << "                 default: no shift, S = 0." << endl;
      cout << "  --invert:      invert mode (look for 1./lambda instead of lambda)." << endl;
      cout << "                 default: no invert" << endl;
      cout << "  --tol T:       tolerance T." << endl;
      cout << "                 default: 1.e-06" << endl;
      cout << "  --maxIt M:     maximum iterations M." << endl;
      cout << "                 default: 100" << endl;
      cout << "  --slv S:       solver (BiCG, CG, LU)" << endl;
      cout << "                   BiCG: iterative method, any matrices" << endl;
      cout << "                     CG: iterative method, sym matrices only" << endl;
      cout << "                     LU:    direct method, any matrices" << endl;
      cout << "                     QR:    direct method, any matrices" << endl;
      cout << "                 default: BiCG" << endl;
      cout << "  --slvTol T:    tolerance T (if iterative method)." << endl;
      cout << "                 default: 1.e-06" << endl;
      cout << "  --slvMaxIt M:  maximum iterations M (if iterative method)." << endl;
      cout << "                 default: 100" << endl;
      cout << "  --noCheck:     check arpack eigen values/vectors." << endl;
      cout << "                 default: check" << endl;
      cout << "  --verbose V:   verbosity level (up to 3)." << endl;
      cout << "                 default: 0" << endl;
      cout << "  --debug D:     debug level (up to 3)." << endl;
      cout << "                 default: 0" << endl;
      cout << "  --restart:     restart from previous run (which had produced resid.out and v.out)." << endl;
      cout << "                 default: false" << endl;
      if (rc == 0) exit(0);
      return rc;
    };

    friend ostream & operator<< (ostream & ostr, options const & opt);

    string fileA;
    string fileB;
    a_int nbEV;
    a_int nbCV;
    bool stdPb; // Standard or generalized (= not standard).
    bool symPb;
    string mag; // Magnitude <=> "which" arpack parameter.
    bool shiftReal, shiftImag;
    double sigmaReal, sigmaImag; // Eigen value translation: look for lambda+sigma instead of lambda.
    bool invert; // Eigen value invertion: look for 1./lambda instead of lambda.
    double tol;
    int maxIt;
    string slv;
    double slvTol;
    int slvMaxIt;
    bool check;
    int verbose;
    int debug;
    bool restart;
};

ostream & operator<< (ostream & ostr, options const & opt) {
  ostr << "OPT: A " << opt.fileA << ", B " << opt.fileB;
  ostr << ", nbEV " << opt.nbEV << ", nbCV " << opt.nbCV << ", stdPb " << (opt.stdPb ? "yes" : "no");
  ostr << ", symPb " << (opt.symPb ? "yes" : "no") << ", mag " << opt.mag << endl;
  ostr << "OPT: shiftReal " << (opt.shiftReal ? "yes" : "no") << ", sigmaReal " << opt.sigmaReal;
  ostr << ", shiftImag " << (opt.shiftImag ? "yes" : "no") << ", sigmaImag " << opt.sigmaImag;
  ostr << ", invert " << (opt.invert ? "yes" : "no") << ", tol " << opt.tol << ", maxIt " << opt.maxIt << endl;
  ostr << "OPT: slv " << opt.slv << ", slvTol " << opt.slvTol << ", slvMaxIt " << opt.slvMaxIt;
  ostr << ", check " << (opt.check ? "yes" : "no") << ", verbose " << opt.verbose << ", debug " << opt.debug;
  ostr << ", restart " << (opt.restart ? "yes" : "no") << endl;
  return ostr;
}

int readMatrixMarket(string const & fileName, EigMatR & M, int const & verbose, string const & msg) {
  ifstream inp(fileName);
  if (!inp) {cerr << "Error: can not open " << fileName << endl; return 1;}

  a_uint l = 0, n = 0, m = 0, nnz = 0;
  vector<a_uint> i, j;
  vector<double> Mij;
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
      double Mkl = 0.;
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

  // Create matrix from file.

  M = EigMatR(n, m); // Set matrice dimensions.
  vector<EigCooR> triplets;
  triplets.reserve(nnz);
  for (size_t k = 0; k < nnz; k++) triplets.emplace_back(i[k], j[k], Mij[k]);
  M.setFromTriplets(triplets.begin(), triplets.end()); // Set all (i, j, Mij).

  if (verbose == 3) {
    cout << endl << msg << endl;
    cout << endl << M << endl;
  }

  return 0;
}

class arpackEV { // Arpack eigen values / vectors.
  public:
    vector<complex<double>> val; // Eigen values.
    vector<EigVecC> vec; // Eigen vectors.
    int nbIt;
    double rciTime;
};

template<typename SLV>
int arpackSolve(options const & opt, int const & mode,
                EigMatR const & A, EigMatR const & B, SLV & solver, arpackEV & out) {
  // Arpack set up.

  // Note: all in/out parameters (all but work*) passed to d[sn][ae]upd are set to 0. before use.
  // d[sn][ae]upd uses dgetv0 to generate a random starting vector (when info is initialized to 0).
  // dgetv0 rely on resid/v: resid/v should be initialized to 0.0 to avoid "bad" starting random vectors.

  char const * which = opt.mag.c_str();
  a_int ido = 0; // First call to arpack.
  char const * iMat = "I";
  char const * gMat = "G";
  char const * bMat = (mode == 1) ? iMat : gMat;
  a_int nbDim = A.rows();
  double * resid = new double[nbDim]; for (a_int n = 0; n < nbDim; n++) resid[n] = 0.; // Avoid "bad" starting vector.
  if (opt.restart) {
    ifstream rfs("resid.out");
    if (rfs.is_open()) {
      for (a_int n = 0; n < nbDim; n++) rfs >> resid[n];
      if (opt.verbose >= 2) {
        cout << endl;
        cout << "resid:" << endl;
        for (a_int n = 0; n < nbDim; n++) cout << resid[n] << endl;
        cout << endl;
      }
    }
  }
  a_int ldv = nbDim;
  double * v = new double[ldv*opt.nbCV]; for (a_int n = 0; n < ldv*opt.nbCV; n++) v[n] = 0.; // Avoid "bad" starting vector.
  if (opt.restart) {
    ifstream vfs("v.out");
    if (vfs.is_open()) {
      a_int nbCV = 0; vfs >> nbCV; if (opt.nbCV < nbCV) nbCV = opt.nbCV;
      for (a_int n = 0; n < ldv*nbCV; n++) vfs >> v[n];
      if (opt.verbose >= 2) {
        cout << endl;
        cout << "v:" << endl;
        for (a_int n = 0; n < ldv*nbCV; n++) cout << v[n] << endl;
        cout << endl;
      }
    }
  }
  a_int iparam[11];
  iparam[0] = 1; // Use exact shifts (=> we'll never have ido == 3).
  iparam[2] = opt.maxIt; // Maximum number of iterations.
  iparam[3] = 1; // Block size.
  iparam[4] = 0; // Number of ev found by arpack.
  if (mode == 1) {
    iparam[6] = mode;
  }
  else if (mode == 2 || mode == 3) {
    if (mode == 2) { // Regular mode.
      iparam[6] = mode;
      solver.compute(B);
    }
    else { // Shift invert mode.
      iparam[6] = mode;
      if (!opt.shiftImag) {
        if (fabs(opt.sigmaReal) < numeric_limits<double>::epsilon()) solver.compute(A);
        else solver.compute(A - opt.sigmaReal * B); // Only real shift.
      }
      else {
        complex<double> sigma(opt.sigmaReal, opt.sigmaImag);
        EigMatC S = A.cast<complex<double>>() - sigma * B.cast<complex<double>>();
        solver.compute(S.real()); // S: shifted matrix is real.
      }
    }
    if(solver.info() != Eigen::Success) {cerr << "Error: decomposition KO - check A and/or B are invertible" << endl; return 1;}
  }
  else {cerr << "Error: arpack mode must be 1, 2 or 3 - KO" << endl; return 1;}
  a_int ipntr[14];
  double * workd = new double[3*nbDim];
  a_int lworkl = opt.symPb ? opt.nbCV*opt.nbCV + 8*opt.nbCV : 3*opt.nbCV*opt.nbCV + 6*opt.nbCV;
  lworkl++; // The documentation says "LWORKL must be at least ..."
  double * workl = new double[lworkl];
  a_int info = 0; // Use random initial residual vector.
  if (opt.restart) info = 1;

  // Arpack solve.

  do {
    // Call arpack.

    if (opt.symPb) {
      dsaupd_c(&ido, bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
      if (info ==  1) cerr << "Error: dsaupd - KO: maximum number of iterations taken. Increase --maxIt..." << endl;
      if (info ==  2) cerr << "Error: dsaupd - KO: no shifts could be applied. Increase --nbCV..." << endl;
      if (info == -9) cerr << "Error: dsaupd - KO: starting vector is zero. Retry: play with shift..." << endl;
      if (info < 0) {cerr << "Error: dsaupd - KO with info " << info << ", nbIt " << iparam[2] << endl; return 1;}
    }
    else {
      dnaupd_c(&ido, bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
      if (info ==  1) cerr << "Error: dnaupd - KO: maximum number of iterations taken. Increase --maxIt..." << endl;
      if (info ==  2) cerr << "Error: dnaupd - KO: no shifts could be applied. Increase --nbCV..." << endl;
      if (info == -9) cerr << "Error: dnaupd - KO: starting vector is zero. Retry: play with shift..." << endl;
      if (info < 0) {cerr << "Error: dnaupd - KO with info " << info << ", nbIt " << iparam[2] << endl; return 1;}
    }

    // Reverse Communication Interface: perform actions according to arpack.

    auto start = chrono::high_resolution_clock::now();

    a_int xIdx = ipntr[0] - 1; // 0-based (Fortran is 1-based).
    a_int yIdx = ipntr[1] - 1; // 0-based (Fortran is 1-based).

    EigMpVR X(workd + xIdx, nbDim); // Arpack provides X.
    EigMpVR Y(workd + yIdx, nbDim); // Arpack provides Y.

    if (ido == -1) {
      if (iparam[6] == 1) {
        Y = A * X;
      }
      else if (iparam[6] == 2) {
        Y = A * X;
        auto YY = Y; // Use copy of Y (not Y) for solve (avoid potential memory overwrite as Y is both in/out).
        Y = solver.solve(YY); // Y = B^-1 * A * X.
        if(solver.info() != Eigen::Success) {cerr << "Error: solve KO - increase --slvMaxIt and/or relax --slvTol, or, change --slv" << endl; return 1;}
      }
      else if (iparam[6] == 3) {
        auto Z = B * X;      // Z = B * X.
        Y = solver.solve(Z); // Y = (A - sigma * B)^-1 * B * X.
        if(solver.info() != Eigen::Success) {cerr << "Error: solve KO - increase --slvMaxIt and/or relax --slvTol, or, change --slv" << endl; return 1;}
      }
    }
    else if (ido == 1) {
      if (iparam[6] == 1) {
        Y = A * X;
      }
      else if (iparam[6] == 2) {
        Y = A * X;
        if (opt.symPb) X = Y; // Remark 5 in dsaupd documentation.
        auto YY = Y; // Use copy of Y (not Y) for solve (avoid potential memory overwrite as Y is both in/out).
        Y = solver.solve(YY);  // Y = B^-1 * A * X.
        if(solver.info() != Eigen::Success) {cerr << "Error: solve KO - increase --slvMaxIt and/or relax --slvTol, or, change --slv" << endl; return 1;}
      }
      else if (iparam[6] == 3) {
        a_int zIdx = ipntr[2] - 1;      // 0-based (Fortran is 1-based).
        EigMpVR Z(workd + zIdx, nbDim); // Arpack provides Z.
        Y = solver.solve(Z);            // Y = (A - sigma * B)^-1 * B * X.
        if(solver.info() != Eigen::Success) {cerr << "Error: solve KO - increase --slvMaxIt and/or relax --slvTol, or, change --slv" << endl; return 1;}
      }
    }
    else if (ido == 2) {
      if      (iparam[6] == 1) Y =     X; // Y = I * X.
      else if (iparam[6] == 2) Y = B * X; // Y = B * X.
      else if (iparam[6] == 3) Y = B * X; // Y = B * X.
    }
    else if (ido != 99) {cerr << "Error: unexpected ido " << ido << " - KO" << endl; return 1;}

    auto stop = chrono::high_resolution_clock::now();
    out.rciTime += chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;

  } while (ido != 99);

  // Get arpack results (computed eigen values and vectors).

  out.nbIt = iparam[2]; // Actual number of iterations.
  bool rvec = true;
  char const * howmny = "A";
  a_int * select = new a_int[opt.nbCV]; for (a_int n = 0; n < opt.nbCV; n++) select[n] = 1;
  a_int const nbZ = nbDim*(opt.nbEV+1); // Caution: opt.nbEV+1 for dneupd.
  double * z = new double[nbZ]; for (a_int n = 0; n < nbZ; n++) z[n] = 0.;
  a_int ldz = nbDim;
  if (opt.symPb) {
    double * d = new double[opt.nbEV]; for (a_int k = 0; k < opt.nbEV; k++) d[k] = 0.;

    dseupd_c(rvec, howmny, select, d, z, ldz, opt.sigmaReal,
             bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
    if (info == -14) cerr << "Error: dseupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
    if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dseupd - KO with info " << info << endl; return 1;}

    // Arpack compute the whole spectrum.

    a_int nbConv = iparam[4];
    out.val.reserve(nbConv);
    for (a_int i = 0; d && i < nbConv; i++) {
      complex<double> lambda(d[i], 0.);
      out.val.push_back(lambda);
      if (out.val.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).
    }

    out.vec.reserve(nbConv);
    for (a_int i = 0; z && i < nbConv; i++) {
      EigVecR V = EigMpVR(z + i*nbDim, nbDim);
      out.vec.push_back(V.cast<complex<double>>());
      if (out.vec.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).
    }

    if (d) {delete [] d; d = NULL;}
  }
  else {
    double * dr = new double[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) dr[k] = 0.;
    double * di = new double[opt.nbEV+1]; for (a_int k = 0; k < opt.nbEV+1; k++) di[k] = 0.;
    double * workev = new double[3*opt.nbCV];

    dneupd_c(rvec, howmny, select, dr, di, z, ldz, opt.sigmaReal, opt.sigmaImag, workev,
             bMat, nbDim, which, opt.nbEV, opt.tol, resid, opt.nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
    if (info == -14) cerr << "Error: dneupd - KO: dnaupd did not find any eigenvalues to sufficient accuracy" << endl;
    if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dneupd - KO with info " << info << endl; return 1;}

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

      EigVecR Vr = EigMpVR(z + (2*i+0)*nbDim, nbDim); // Real part.
      EigVecR Vi = EigMpVR(z + (2*i+1)*nbDim, nbDim); // Imaginary part.
      complex<double> imag(0., 1.);
      EigVecC V = Vr.cast<complex<double>>() + imag * Vi.cast<complex<double>>();
      out.vec.push_back(V);
      if (out.vec.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).

      // Deduce second half of the spectrum.

      V = Vr.cast<complex<double>>() - imag * Vi.cast<complex<double>>();
      out.vec.push_back(V);
      if (out.vec.size() == (size_t) opt.nbEV) break; // If more converged than requested, likely not accurate (check KO).
    }

    if (workev) {delete [] workev; workev = NULL;}
    if (dr)     {delete [] dr;     dr     = NULL;}
    if (di)     {delete [] di;     di     = NULL;}
  }

  ofstream rfs("resid.out"); for (a_int n = 0; n < nbDim; n++) rfs << resid[n] << endl;
  ofstream vfs("v.out"); vfs << opt.nbCV << endl; for (a_int n = 0; n < ldv*opt.nbCV; n++) vfs << v[n] << endl;

  // Clean.

  if (z)      {delete [] z;      z      = NULL;}
  if (select) {delete [] select; select = NULL;}
  if (workl)  {delete [] workl;  workl  = NULL;}
  if (workd)  {delete [] workd;  workd  = NULL;}
  if (v)      {delete [] v;      v      = NULL;}
  if (resid)  {delete [] resid;  resid  = NULL;}

  return 0;
}

int checkArpackEigVec(options const & opt, EigMatR & A, EigMatR const & B, arpackEV const & out) {
  // Check eigen vectors.

  if (opt.check && out.vec.size() == 0) {
    cerr << "Error: no eigen value / vector found" << endl;
    return 1;
  }

  for (size_t i = 0; i < out.vec.size(); i++) {
    EigVecC V = out.vec[i];
    complex<double> lambda = out.val[i];
    if (opt.verbose >= 1) {
      cout << endl;
      cout << "eigen value " << setw(3) << i << ": " << lambda << endl;
      if (opt.verbose >= 2) {
        cout << endl;
        cout << "eigen vector " << setw(3) << i << " (norm " << V.norm() << "): " << endl;
        cout << endl << V << endl;
      }
    }

    if (opt.check) {
      EigVecC left = A.cast<complex<double>>() * V;
      EigVecC right = opt.stdPb ? V : B.cast<complex<double>>() * V;
      right *= lambda;
      EigVecC diff = left - right;
      if (diff.norm() > sqrt(opt.tol)) {
        cerr << endl << "Error: bad eigen vector " << setw(3) << i << " (norm " << V.norm() << "):" << endl;
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
          cout << endl << "eigen value/vector " << setw(3) << i << ": check OK";
          cout << ", diff (norm " << diff.norm() << ", sqrt(tol) " << sqrt(opt.tol) << ")" << endl;
        }
      }
    }
  }

  return 0;
}

template<typename SLV>
int arpackSolve(options const & opt, EigMatR & A, EigMatR const & B,
                SLV & solver, arpackEV & out) {
  // If needed, transform the initial problem into a new one that arpack can handle.

  auto eps = numeric_limits<double>::epsilon();
  bool shiftReal = (opt.shiftReal && fabs(opt.sigmaReal) > eps) ? true : false;
  bool shiftImag = (opt.shiftImag && fabs(opt.sigmaImag) > eps) ? true : false;

  bool backTransform = false;
  int mode = 0;
  if (opt.stdPb) {
    mode = 1;
    if (shiftReal && !shiftImag) {
      EigMatR I(A.rows(), A.cols());
      I.setIdentity();
      A -= opt.sigmaReal*I;
      backTransform = true;
    }
  }
  else {
    mode = 2;
    if (shiftReal || shiftImag) mode = 3;
  }

  // Solve the problem.

  if (opt.verbose >= 1) {
    cout << endl;
    cout << "ARP: mode " << mode;
    cout << ", nbDim " << A.rows();
    cout << ", backTransform " << (backTransform ? "yes" : "no") << endl;
  }

  int rc = arpackSolve<SLV>(opt, mode, A, B, solver, out);
  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

  if (opt.verbose >= 1) {
    cout << endl;
    cout << "ARP: nbEV found " << out.val.size();
    cout << ", nbIt " << out.nbIt << endl;
  }

  // If needed, transform back the arpack problem into the initial problem.

  if (backTransform) {
    for (size_t i = 0; i < out.val.size(); i++) out.val[i] += opt.sigmaReal;
    EigMatR I(A.rows(), A.cols());
    I.setIdentity();
    A += opt.sigmaReal*I; // For later checks.
  }

  // Check.

  return checkArpackEigVec(opt, A, B, out);
}

template<typename SLV>
int arpackSolve(options & opt, SLV & solver) {
  // Read A.

  EigMatR A;
  int rc = readMatrixMarket(opt.fileA, A, opt.verbose, "A:");
  if (rc != 0) {cerr << "Error: read A KO" << endl; return rc;}

  // Read B.

  EigMatR B;
  if (!opt.stdPb) {
    rc = readMatrixMarket(opt.fileB, B, opt.verbose, "B:");
    if (rc != 0) {cerr << "Error: read B KO" << endl; return rc;}
  }

  // Check A-B compatibility.

  if (!opt.stdPb) {
    if (A.rows() != B.rows()) {cerr << "Error: A.rows() != B.rows()" << endl; return rc;}
    if (A.cols() != B.cols()) {cerr << "Error: A.cols() != B.cols()" << endl; return rc;}
  }
  if (opt.nbCV > A.cols()) opt.nbCV = A.cols(); // Cut-off.

  // Arpack solve.

  arpackEV out;
  out.rciTime = 0.;
  auto start = chrono::high_resolution_clock::now();
  rc = arpackSolve<SLV>(opt, A, B, solver, out);
  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}
  auto stop = chrono::high_resolution_clock::now();
  double fullTime = chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;
  cout << endl;
  cout << "OUT: nb EV found " << out.val.size() << ", nb iterations " << out.nbIt << endl;
  cout << "OUT: full time " << fullTime << " s, RCI time " << out.rciTime << " s" << endl;

  return 0;
}

int main(int argc, char ** argv) {
  // Check for options.

  options opt;
  int rc = opt.readCmdLine(argc, argv);
  if (rc != 0) {cerr << "Error: read cmd line KO" << endl; return rc;}
  cout << opt; // Print options.

  // Solve with arpack.

  if (opt.slv == "BiCG") {
    EigBiCG solver;
    solver.setTolerance(opt.slvTol);
    solver.setMaxIterations(opt.slvMaxIt);
    rc = arpackSolve<EigBiCG>(opt, solver);
  }
  else if (opt.slv == "CG") {
    EigCG solver;
    solver.setTolerance(opt.slvTol);
    solver.setMaxIterations(opt.slvMaxIt);
    rc = arpackSolve<EigCG>(opt, solver);
  }
  else if (opt.slv == "LU") {
    EigSLU solver;
    rc = arpackSolve<EigSLU>(opt, solver);
  }
  else if (opt.slv == "QR") {
    EigSQR solver;
    rc = arpackSolve<EigSQR>(opt, solver);
  }
  else {cerr << "Error: unknown solver - KO" << endl; return 1;}
  if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

  return 0;
}

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=2 ts=2 et smartindent :*/
