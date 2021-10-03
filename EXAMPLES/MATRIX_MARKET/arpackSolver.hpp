#ifndef __ARPACKSOLVER_HPP__
#define __ARPACKSOLVER_HPP__

#include "arpack.h"

#include <iostream>
#include <string>
#include <sstream> // stringstream.
#include <fstream> // [io]fstream.
#include <chrono>
#include <iomanip> // setw.
#include <cassert>
#include <vector>
#include <type_traits> // is_same.

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

// Sparse matrix related types.

typedef Eigen::SparseMatrix<         float>  EigSMxS; // Real.
typedef Eigen::SparseMatrix<        double>  EigSMxD; // Real.
typedef Eigen::SparseMatrix<complex< float>> EigSMxC; // Complex.
typedef Eigen::SparseMatrix<complex<double>> EigSMxZ; // Complex.

// Iterative solvers for sparse matrices.

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

// Direct solvers for sparse matrices.

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

// Direct solvers for dense matrices.

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

// Definition of arpackSolver class.

typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, 1> EigVecZ;

typedef vector<complex<double>> StdVecZ;
typedef vector<EigVecZ> StdVecEVZ;

// RC: Real or Complex.
// FD: Float or Double.
// EM: Eigen Matrix (sparse or dense).
// SLV: Solver.
template<typename RC, typename FD, typename EM, typename SLV>
class arpackSolver {
  // Nested typedef.

  typedef Eigen::Map<Eigen::Matrix<RC, Eigen::Dynamic, 1>> EV;

  // Public methods.

  public:

    arpackSolver() {
      symPb = true;
      nbEV = 1;
      nbCV = 2*nbEV+1;
      tol = 1.e-6;
      sigmaReal = sigmaImag = 0.;
      dumpToFile = false;
      restartFromFile = false;
      mag = "LM";
      maxIt = 100;
      schur = false;
      verbose = 0;

      stdPb = true;
      mode = 1;
      nbIt = 0;
      imsTime = 0.;
      rciTime = 0.;

      nbDim = 0;
      resid = nullptr;
      v = nullptr;
    };

    ~arpackSolver() {
      if (v)     {delete [] v;     v      = nullptr;}
      if (resid) {delete [] resid; resid  = nullptr;}
    };

    int createMatrix(string const & fileName, Eigen::SparseMatrix<RC> & M) {
      // Read matrix from file.

      a_uint n = 0, m = 0;
      vector<a_uint> i, j;
      vector<RC> Mij;
      int rc = readMatrixMarket(fileName, n, m, i, j, Mij);
      if (rc != 0) {cerr << "Error: read matrix market file KO" << endl; return rc;}

      // Create matrix from file.

      M = Eigen::SparseMatrix<RC>(n, m); // Set matrice dimensions.
      vector<Eigen::Triplet<RC>> triplets;
      a_uint nnz = Mij.size();
      triplets.reserve(nnz);
      for (size_t k = 0; k < nnz; k++) triplets.emplace_back(i[k], j[k], Mij[k]);
      M.setFromTriplets(triplets.begin(), triplets.end()); // Set all (i, j, Mij).

      return 0;
    };

    int createMatrix(string const & fileName, Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic> & M) {
      // Read matrix from file.

      a_uint n = 0, m = 0;
      vector<a_uint> i, j;
      vector<RC> Mij;
      int rc = readMatrixMarket(fileName, n, m, i, j, Mij);
      if (rc != 0) {cerr << "Error: read matrix market file KO" << endl; return rc;}

      // Create matrix from file.

      M = Eigen::Matrix<RC, Eigen::Dynamic, Eigen::Dynamic>(n, m); // Set matrice dimensions.
      M.setZero(n, m); // Avoid spurious/random values which may break solves (LU, QR, ...).
      a_uint nnz = Mij.size();
      for (size_t k = 0; k < nnz; k++) M(i[k], j[k]) = Mij[k];

      return 0;
    };

    void dumpParameters() {
      if (verbose >= 1) {
        cout << endl << "arpackSolver:" << endl;
        cout << endl << "symPb: " << symPb << endl;
        cout << endl << "nbEV: " << nbEV << endl;
        cout << endl << "nbCV: " << nbCV << endl;
        cout << endl << "tol: " << tol << endl;
        cout << endl << "sigmaReal: " << sigmaReal << endl;
        cout << endl << "sigmaImag: " << sigmaImag << endl;
        cout << endl << "dumpToFile: " << dumpToFile << endl;
        cout << endl << "restartFromFile: " << restartFromFile << endl;
        cout << endl << "mag: " << mag << endl;
        cout << endl << "maxIt: " << maxIt << endl;
        cout << endl << "schur: " << schur << endl;
      }
    };

    int solve(EM & A, EM const * B = nullptr) {
      stdPb = !B ? true : false;

      dumpParameters();
      if (verbose == 3) {
        cout << endl << "arpackSolver:" << endl;
        cout << endl << "A:" << endl;
        cout << endl <<  A   << endl;
        if (B) {
          cout << endl << "B:" << endl;
          cout << endl << *B   << endl;
        }
      }
      nbDim = A.rows();

      // If needed, transform the initial problem into a new one that arpack can handle.

      auto eps = numeric_limits<double>::epsilon();
      bool shiftReal = (fabs(sigmaReal) > eps) ? true : false;
      bool shiftImag = (fabs(sigmaImag) > eps) ? true : false;
      bool backTransform = false;
      mode = 0;
      if (stdPb) {
        mode = 1;
        if (shiftReal && !shiftImag) {
          EM I(A.rows(), A.cols());
          I.setIdentity();
          RC sigma; makeSigma(sigma);
          A -= sigma*I;
          backTransform = true;
        }
      }
      else {
        mode = 2;
        if (shiftReal || shiftImag) mode = 3;
      }
      if (verbose >= 1) {
        cout << endl << "arpackSolver:" << endl;
        cout << endl << "mode " << mode << ", backTransform " << (backTransform ? "yes" : "no") << endl;
      }

      // Solve with arpack.
      // Note: when initializing Eigen solvers, API differ depending on solvers.

      SLV solver;
      int rc = initSolver(solver);
      if (rc != 0) {cerr << "Error: initialize solver KO" << endl; return rc;}
      rc = solve(A, B, solver);
      if (rc != 0) {cerr << "Error: arpack solve KO" << endl; return rc;}

      // If needed, transform back the arpack problem into the initial problem.

      if (backTransform) {
        EM I(A.rows(), A.cols());
        I.setIdentity();
        RC sigma; makeSigma(sigma);
        A += sigma*I;
        for (size_t i = 0; i < val.size(); i++) val[i] += sigma;
      }

      return 0;
    };

    int checkEigVec(EM const & A, EM const * B = nullptr, double const * diffTol = nullptr) {
      stdPb = !B ? true : false;
      double dTol = !diffTol ? sqrt(tol) : *diffTol;

      // Check eigen vectors.

      string rs = schur ? "Schur" : "Ritz";

      if (vec.size() == 0) {
        cerr << "Error: no " << rs << " value / vector found" << endl;
        return 1;
      }

      for (size_t i = 0; i < vec.size(); i++) {
        EigVecZ V = vec[i];
        complex<double> lambda = val[i];
        if (verbose >= 1) {
          cout << endl << "arpackSolver:" << endl;
          cout << endl << rs << " value " << setw(3) << i << ": " << lambda << endl;
          if (verbose >= 2) {
            cout << endl << rs << " vector " << setw(3) << i << " (norm " << V.norm() << "): " << endl;
            cout << endl << V << endl;
          }
        }

        EigVecZ left = A.template cast<complex<double>>() * V;
        EigVecZ right = stdPb ? V : B->template cast<complex<double>>() * V;
        right *= lambda;
        EigVecZ diff = left - right;
        if (diff.norm() > dTol) {
          cerr << endl << "Error: bad vector " << setw(3) << i << " (norm " << V.norm() << "):" << endl;
          cerr << endl << V << endl;
          cerr << endl << "Error: left side (A*V - norm " << left.norm() << "):" << endl;
          cerr << endl << left << endl;
          cerr << endl << "Error: right side (lambda*" << (stdPb ? "" : "B*") << "V - norm " << right.norm() << "):" << endl;
          cerr << endl << right << endl;
          cerr << endl << "Error: diff (norm " << diff.norm() << ", tol " << dTol << "):" << endl;
          cerr << endl << diff << endl;
          return 1;
        }
        else {
          if (verbose >= 1) {
            cout << endl << "arpackSolver:" << endl;
            cout << endl << rs << " value/vector " << setw(3) << i << ": check OK";
            cout << ", diff (norm " << diff.norm() << ", tol " << dTol << ")" << endl;
          }
        }
      }

      return 0;
    };

  // Private methods.

  private:

    void makeZero(         float  & zero) {zero = 0.f;};

    void makeZero(        double  & zero) {zero = 0.;};

    void makeZero(complex< float> & zero) {zero = complex<double>(0.f, 0.f);};

    void makeZero(complex<double> & zero) {zero = complex<double>(0., 0.);};

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
    };

    void aupd(a_int * ido, char const * bMat, a_int nbDim, char const * which, float * resid, float * v,
              a_int ldv, a_int * iparam, a_int * ipntr, float * workd, float * workl, a_int lworkl, float * & rwork,
              a_int * info) {
      assert(rwork == nullptr);
      if (symPb) {
        ssaupd_c(ido, bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
      }
      else {
        snaupd_c(ido, bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
      }
    };

    void aupd(a_int * ido, char const * bMat, a_int nbDim, char const * which, double * resid, double * v,
              a_int ldv, a_int * iparam, a_int * ipntr, double * workd, double * workl, a_int lworkl, double * & rwork,
              a_int * info) {
      assert(rwork == nullptr);
      if (symPb) {
        dsaupd_c(ido, bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
      }
      else {
        dnaupd_c(ido, bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, info);
      }
    };

    void aupd(a_int * ido, char const * bMat, a_int nbDim, char const * which, complex<float> * resid, complex<float> * v,
              a_int ldv, a_int * iparam, a_int * ipntr, complex<float> * workd, complex<float> * workl, a_int lworkl, float * & rwork,
              a_int * info) {
      if (!rwork) rwork = new float[nbCV];
      cnaupd_c(ido, bMat, nbDim, which, nbEV, tol, reinterpret_cast<_Complex float*>(resid), nbCV,
               reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr, reinterpret_cast<_Complex float*>(workd),
               reinterpret_cast<_Complex float*>(workl), lworkl, rwork, info);
    };

    void aupd(a_int * ido, char const * bMat, a_int nbDim, char const * which, complex<double> * resid, complex<double> * v,
              a_int ldv, a_int * iparam, a_int * ipntr, complex<double> * workd, complex<double> * workl, a_int lworkl, double * & rwork,
              a_int * info) {
      if (!rwork) rwork = new double[nbCV];
      znaupd_c(ido, bMat, nbDim, which, nbEV, tol, reinterpret_cast<_Complex double*>(resid), nbCV,
               reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr, reinterpret_cast<_Complex double*>(workd),
               reinterpret_cast<_Complex double*>(workl), lworkl, rwork, info);
    };

    void spectrum(RC * d, RC * z, a_int nbDim, a_int * iparam) {
      // Arpack compute the whole spectrum.

      a_int nbConv = iparam[4];
      val.reserve(nbConv);
      for (a_int i = 0; d && i < nbConv; i++) {
        complex<double> lambda(d[i]);
        val.push_back(lambda);
        if (val.size() == (size_t) nbEV) break; // If more converged than requested, likely not accurate (check KO).
      }

      vec.reserve(nbConv);
      for (a_int i = 0; z && i < nbConv; i++) {
        EV V = EV(z + i*nbDim, nbDim);
        vec.push_back(V.template cast<complex<double>>());
        if (vec.size() == (size_t) nbEV) break; // If more converged than requested, likely not accurate (check KO).
      }
    };

    void halfSpectrum(RC * dr, RC * di, RC * z, a_int nbDim, a_int * iparam) {
      // Arpack compute only half of the spectrum.

      a_int nbConv = iparam[4];
      val.reserve(nbConv);
      for (a_int i = 0; dr && di && i <= nbConv/2; i++) { // Scan first half of the spectrum.
        // Get first half of the spectrum.

        complex<double> lambda(dr[i], di[i]);
        val.push_back(lambda);
        if (val.size() == (size_t) nbEV) break; // If more converged than requested, likely not accurate (check KO).

        // Deduce second half of the spectrum.

        val.push_back(complex<double>(lambda.real(), -1.*lambda.imag()));
        if (val.size() == (size_t) nbEV) break; // If more converged than requested, likely not accurate (check KO).
      }

      vec.reserve(nbConv);
      for (a_int i = 0; z && i <= nbConv/2; i++) { // Scan half spectrum.
        // Get first half of the spectrum.

        EV Vr = EV(z + (2*i+0)*nbDim, nbDim); // Real part.
        EV Vi = EV(z + (2*i+1)*nbDim, nbDim); // Imaginary part.
        complex<double> imag(0., 1.);
        EigVecZ V = Vr.template cast<complex<double>>() + imag * Vi.template cast<complex<double>>();
        vec.push_back(V);
        if (vec.size() == (size_t) nbEV) break; // If more converged than requested, likely not accurate (check KO).

        // Deduce second half of the spectrum.

        V = Vr.template cast<complex<double>>() - imag * Vi.template cast<complex<double>>();
        vec.push_back(V);
        if (vec.size() == (size_t) nbEV) break; // If more converged than requested, likely not accurate (check KO).
      }
    };

    int eupd(a_int rvec, char const * howmny, a_int const * select, float * z,
             a_int ldz, char const * bMat, a_int nbDim, char const * which, float * resid, float * v,
             a_int ldv, a_int * iparam, a_int * ipntr, float * workd, float * workl, a_int lworkl, float * rwork,
             a_int & info) {
      assert(rwork == nullptr);
      if (symPb) {
        float * d = new float[nbEV]; for (a_int k = 0; k < nbEV; k++) d[k] = 0.;

        sseupd_c(rvec, howmny, select, d, z, ldz, sigmaReal,
                 bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
        if (info == -14) cerr << "Error: dseupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
        if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dseupd - KO with info " << info << endl; return 1;}

        spectrum(d, z, nbDim, iparam);

        if (d) {delete [] d; d = nullptr;}
      }
      else {
        float * dr = new float[nbEV+1]; for (a_int k = 0; k < nbEV+1; k++) dr[k] = 0.;
        float * di = new float[nbEV+1]; for (a_int k = 0; k < nbEV+1; k++) di[k] = 0.;
        float * workev = new float[3*nbCV];

        sneupd_c(rvec, howmny, select, dr, di, z, ldz, sigmaReal, sigmaImag, workev,
                 bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
        if (info == -14) cerr << "Error: dneupd - KO: [dz]naupd did not find any eigenvalues to sufficient accuracy" << endl;
        if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dneupd - KO with info " << info << endl; return 1;}

        halfSpectrum(dr, di, z, nbDim, iparam);

        if (workev) {delete [] workev; workev = nullptr;}
        if (dr)     {delete [] dr;     dr     = nullptr;}
        if (di)     {delete [] di;     di     = nullptr;}
      }

      return 0;
    };

    int eupd(a_int rvec, char const * howmny, a_int const * select, double * z,
             a_int ldz, char const * bMat, a_int nbDim, char const * which, double * resid, double * v,
             a_int ldv, a_int * iparam, a_int * ipntr, double * workd, double * workl, a_int lworkl, double * rwork,
             a_int & info) {
      assert(rwork == nullptr);
      if (symPb) {
        double * d = new double[nbEV]; for (a_int k = 0; k < nbEV; k++) d[k] = 0.;

        dseupd_c(rvec, howmny, select, d, z, ldz, sigmaReal,
                 bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
        if (info == -14) cerr << "Error: dseupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
        if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dseupd - KO with info " << info << endl; return 1;}

        spectrum(d, z, nbDim, iparam);

        if (d) {delete [] d; d = nullptr;}
      }
      else {
        double * dr = new double[nbEV+1]; for (a_int k = 0; k < nbEV+1; k++) dr[k] = 0.;
        double * di = new double[nbEV+1]; for (a_int k = 0; k < nbEV+1; k++) di[k] = 0.;
        double * workev = new double[3*nbCV];

        dneupd_c(rvec, howmny, select, dr, di, z, ldz, sigmaReal, sigmaImag, workev,
                 bMat, nbDim, which, nbEV, tol, resid, nbCV, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
        if (info == -14) cerr << "Error: dneupd - KO: [dz]naupd did not find any eigenvalues to sufficient accuracy" << endl;
        if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: dneupd - KO with info " << info << endl; return 1;}

        halfSpectrum(dr, di, z, nbDim, iparam);

        if (workev) {delete [] workev; workev = nullptr;}
        if (dr)     {delete [] dr;     dr     = nullptr;}
        if (di)     {delete [] di;     di     = nullptr;}
      }

      return 0;
    };

    int eupd(a_int rvec, char const * howmny, a_int const * select, complex<float> * z,
             a_int ldz, char const * bMat, a_int nbDim, char const * which, complex<float> * resid, complex<float> * v,
             a_int ldv, a_int * iparam, a_int * ipntr, complex<float> * workd, complex<float> * workl, a_int lworkl, float * rwork,
             a_int & info) {
      complex<float> * d = new complex<float>[nbEV+1]; for (a_int k = 0; k < nbEV+1; k++) d[k] = complex<float>(0., 0.);
      complex<float> * workev = new complex<float>[2*nbCV];

      complex<float> sigma = complex<float>((float) sigmaReal, (float) sigmaImag);
      cneupd_c(rvec, howmny, select, reinterpret_cast<_Complex float*>(d), reinterpret_cast<_Complex float*>(z), ldz,
               reinterpret_cast<_Complex float &>(sigma), reinterpret_cast<_Complex float*>(workev),
               bMat, nbDim, which, nbEV, tol, reinterpret_cast<_Complex float*>(resid), nbCV,
               reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
               reinterpret_cast<_Complex float*>(workd), reinterpret_cast<_Complex float*>(workl), lworkl, rwork, &info);
      if (info == -14) cerr << "Error: zneupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
      if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: zneupd - KO with info " << info << endl; return 1;}

      spectrum(d, z, nbDim, iparam);

      if (workev) {delete [] workev; workev = nullptr;}
      if (d)      {delete [] d;           d = nullptr;}

      return 0;
    };

    int eupd(a_int rvec, char const * howmny, a_int const * select, complex<double> * z,
             a_int ldz, char const * bMat, a_int nbDim, char const * which, complex<double> * resid, complex<double> * v,
             a_int ldv, a_int * iparam, a_int * ipntr, complex<double> * workd, complex<double> * workl, a_int lworkl, double * rwork,
             a_int & info) {
      complex<double> * d = new complex<double>[nbEV+1]; for (a_int k = 0; k < nbEV+1; k++) d[k] = complex<double>(0., 0.);
      complex<double> * workev = new complex<double>[2*nbCV];

      complex<double> sigma = complex<double>(sigmaReal, sigmaImag);
      zneupd_c(rvec, howmny, select, reinterpret_cast<_Complex double*>(d), reinterpret_cast<_Complex double*>(z), ldz,
               reinterpret_cast<_Complex double &>(sigma), reinterpret_cast<_Complex double*>(workev),
               bMat, nbDim, which, nbEV, tol, reinterpret_cast<_Complex double*>(resid), nbCV,
               reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
               reinterpret_cast<_Complex double*>(workd), reinterpret_cast<_Complex double*>(workl), lworkl, rwork, &info);
      if (info == -14) cerr << "Error: zneupd - KO: dsaupd did not find any eigenvalues to sufficient accuracy" << endl;
      if (info < 0 && info != -14 /*-14: don't break*/) {cerr << "Error: zneupd - KO with info " << info << endl; return 1;}

      spectrum(d, z, nbDim, iparam);

      if (workev) {delete [] workev; workev = nullptr;}
      if (d)      {delete [] d;           d = nullptr;}

      return 0;
    };

    int setMode(int const mode, EM const & A, EM const * B, SLV & solver) {
      int rc = 1;

      if (mode == 1) {
        rc = 0;
      }
      else if (mode == 2 || mode == 3) {
        if (!stdPb && !B) {cerr << "Error: generalized problem without B" << endl; return 1;}

        if (mode == 2) { // Invert mode.
          solver.compute(*B);
        }
        else { // Shift invert mode.
          RC sigma; makeSigma(sigma);
          auto S = A - sigma*(*B);
          solver.compute(S);
        }

        rc = 0;
      }
      else {cerr << "Error: arpack mode must be 1, 2 or 3 - KO" << endl; rc = 1;}

      return rc;
    };

    void restartSolve(string const & fileName,
                      a_int const & nbDim, RC * rv, bool readNbCV) {
      if (restartFromFile) {
        ifstream ifs(fileName.c_str());
        if (ifs.is_open()) {
          a_int nbCV = 1;
          if (readNbCV) {
            ifs >> nbCV;
            if (nbCV < nbDim) nbCV = nbDim; // Cut-off in case previous run used different nbDim.
          }
          for (a_int n = 0; rv && n < nbDim*nbCV; n++) ifs >> rv[n];
          if (verbose >= 1) {
            cout << endl << "arpackSolver:" << endl;
            cout << endl << fileName << ": restart OK" << endl;
            if (verbose >= 2) {
              for (a_int n = 0; rv && n < nbDim; n++) cout << rv[n] << endl;
            }
          }
        }
      }
    };

    int initPointerSize(a_int & iparamSz, a_int & ipntrSz, string const & aeupd) {
      iparamSz = 11; ipntrSz = 14;
      if (aeupd == "aupd") {
        if (is_same<RC, double>::value && symPb) ipntrSz = 11;
        if (is_same<RC,  float>::value && symPb) ipntrSz = 11;
        return 0;
      }
      if (aeupd == "eupd") {
        if (is_same<RC, double>::value && symPb) {iparamSz = 7; ipntrSz = 11;}
        if (is_same<RC,  float>::value && symPb) {iparamSz = 7; ipntrSz = 11;}
        return 0;
      }
      return 1;
    };

    int solve(EM const & A, EM const * B, SLV & solver) {
      if (!stdPb && !B) {cerr << "Error: generalized problem without B" << endl; return 1;}

      // Arpack set up.

      // Note: some in/out parameters passed to d[sn][ae]upd are set to 0. before use.
      // d[sn][ae]upd uses dgetv0 to generate a random starting vector (when info is initialized to 0).
      // dgetv0 rely on resid/v: resid/v should be initialized to 0.0 to avoid "bad" starting random vectors.

      char const * which = mag.c_str();
      a_int ido = 0; // First call to arpack.
      char const * iMat = "I";
      char const * gMat = "G";
      char const * bMat = (mode == 1) ? iMat : gMat;
      a_int nbDim = A.rows();
      RC zero; makeZero(zero);
      if (!resid) {
        resid = new RC[nbDim];
        for (a_int n = 0; n < nbDim; n++) resid[n] = zero; // Avoid "bad" starting vector.
      };
      restartSolve("arpackSolver.resid.out", nbDim, resid, false);
      a_int ldv = nbDim;
      if (!v) {
        v = new RC[ldv*nbCV];
        for (a_int n = 0; n < ldv*nbCV; n++) v[n] = zero; // Avoid "bad" starting vector.
      };
      restartSolve("arpackSolver.v.out", ldv, v, true);
      a_int iparamSz = 0, ipntrSz = 0;
      int rc = initPointerSize(iparamSz, ipntrSz, "aupd");
      if (rc != 0) {cerr << "Error: bad iparam/ipntr initialization for aupd" << endl; return rc;}
      vector<a_int> iparamAupd(iparamSz, 0);
      iparamAupd[0] = 1; // Use exact shifts (=> we'll never have ido == 3).
      iparamAupd[2] = maxIt; // Maximum number of iterations.
      iparamAupd[3] = 1; // Block size.
      iparamAupd[4] = 0; // Number of ev found by arpack.
      iparamAupd[6] = mode;
      vector<a_int> ipntrAupd(ipntrSz, 0);
      RC * workd = new RC[3*nbDim]; for (a_int n = 0; n < 3*nbDim; n++) workd[n] = zero; // Avoid "bad" X/Y vector.
      a_int lworkl = symPb ? nbCV*nbCV + 8*nbCV : 3*nbCV*nbCV + 6*nbCV;
      RC * workl = new RC[lworkl];
      a_int info = 0; // Use random initial residual vector.
      if (restartFromFile) info = 1;

      // Initialize solver.

      auto start = chrono::high_resolution_clock::now();
      rc = setMode(mode, A, B, solver);
      if (rc != 0) {cerr << "Error: bad arpack mode" << endl; return rc;}
      auto stop = chrono::high_resolution_clock::now();
      imsTime = chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;

      // Arpack solve.

      FD * rwork = nullptr;
      do {
        // Call arpack.

        aupd(&ido, bMat, nbDim, which, resid, v, ldv, iparamAupd.data(), ipntrAupd.data(), workd, workl, lworkl, rwork, &info);
        if (info ==  1) cerr << "Error: [dz][sn]aupd - KO: maximum number of iterations taken. Increase --maxIt..." << endl;
        if (info ==  3) cerr << "Error: [dz][sn]aupd - KO: no shifts could be applied. Increase --nbCV..." << endl;
        if (info == -9) cerr << "Error: [dz][sn]aupd - KO: starting vector is zero. Retry: play with shift..." << endl;
        if (info < 0) {cerr << "Error: [dz][sn]aupd - KO with info " << info << ", nbIt " << iparamAupd[2] << endl; return 1;}

        // Reverse Communication Interface: perform actions according to arpack.

        start = chrono::high_resolution_clock::now();

        a_int xIdx = ipntrAupd[0] - 1; // 0-based (Fortran is 1-based).
        a_int yIdx = ipntrAupd[1] - 1; // 0-based (Fortran is 1-based).

        EV X(workd + xIdx, nbDim); // Arpack provides X.
        EV Y(workd + yIdx, nbDim); // Arpack provides Y.

        if (ido == -1) {
          if (iparamAupd[6] == 1) {
            Y = A * X;
          }
          else if (iparamAupd[6] == 2) {
            Y = A * X;
            auto YY = Y;          // Use copy of Y (not Y) for solve (avoid potential memory overwrite as Y is both in/out).
            Y = solver.solve(YY); // Y = B^-1 * A * X.
          }
          else if (iparamAupd[6] == 3) {
            auto Z = (*B) * X;   // Z = B * X.
            Y = solver.solve(Z); // Y = (A - sigma * B)^-1 * B * X.
          }
        }
        else if (ido == 1) {
          if (iparamAupd[6] == 1) {
            Y = A * X;
          }
          else if (iparamAupd[6] == 2) {
            Y = A * X;
            if (symPb) X = Y;     // Remark 5 in dsaupd documentation.
            auto YY = Y;          // Use copy of Y (not Y) for solve (avoid potential memory overwrite as Y is both in/out).
            Y = solver.solve(YY); // Y = B^-1 * A * X.
          }
          else if (iparamAupd[6] == 3) {
            a_int zIdx = ipntrAupd[2] - 1; // 0-based (Fortran is 1-based).
            EV Z(workd + zIdx, nbDim); // Arpack provides Z.
            Y = solver.solve(Z);       // Y = (A - sigma * B)^-1 * B * X.
          }
        }
        else if (ido == 2) {
          if      (iparamAupd[6] == 1) Y =        X; // Y = I * X.
          else if (iparamAupd[6] == 2) Y = (*B) * X; // Y = B * X.
          else if (iparamAupd[6] == 3) Y = (*B) * X; // Y = B * X.
        }
        else if (ido != 99) {cerr << "Error: unexpected ido " << ido << " - KO" << endl; return 1;}

        stop = chrono::high_resolution_clock::now();
        rciTime += chrono::duration_cast<chrono::milliseconds>(stop - start).count()/1000.;

      } while (ido != 99);

      // Get arpack results (computed eigen values and vectors).

      nbIt = iparamAupd[2]; // Actual number of iterations.
      a_int rvec = 1;
      char const * howmnyA = "A"; // Ritz vectors.
      char const * howmnyP = "P"; // Schur vectors.
      char const * howmny = schur ? howmnyP : howmnyA;
      a_int * select = new a_int[nbCV]; for (a_int n = 0; n < nbCV; n++) select[n] = 1;
      a_int const nbZ = nbDim*(nbEV+1); // Caution: nbEV+1 for dneupd.
      RC * z = new RC[nbZ]; for (a_int n = 0; n < nbZ; n++) z[n] = zero;
      a_int ldz = nbDim;
      rc = initPointerSize(iparamSz, ipntrSz, "eupd");
      if (rc != 0) {cerr << "Error: bad iparam/ipntr initialization for eupd" << endl; return rc;}
      vector<a_int> iparamEupd(iparamSz, 0);
      for (a_int p = 0; p < iparamSz; p++) iparamEupd[p] = iparamAupd[p]; // Initialize eupd parameters with aupd ones.
      vector<a_int> ipntrEupd(ipntrSz, 0);
      for (a_int p = 0; p < ipntrSz; p++) ipntrEupd[p] = ipntrAupd[p]; // Initialize eupd parameters with aupd ones.
      rc = eupd(rvec, howmny, select, z, ldz, bMat, nbDim, which, resid, v, ldv, iparamEupd.data(), ipntrEupd.data(), workd, workl, lworkl, rwork, info);
      if (rc != 0) {cerr << "Error: bad arpack eupd" << endl; return rc;}

      if (dumpToFile) {
        ofstream rfs("arpackSolver.resid.out"); for (a_int n = 0; n < nbDim; n++) rfs << resid[n] << endl;
        ofstream vfs("arpackSolver.v.out"); vfs << nbCV << endl; for (a_int n = 0; n < ldv*nbCV; n++) vfs << v[n] << endl;
      }

      // Clean.

      if (rwork)  {delete [] rwork;  rwork  = nullptr;}
      if (z)      {delete [] z;      z      = nullptr;}
      if (select) {delete [] select; select = nullptr;}
      if (workl)  {delete [] workl;  workl  = nullptr;}
      if (workd)  {delete [] workd;  workd  = nullptr;}

      return 0;
    };

    void makeSigma(         float  & sigma) {sigma = (float) sigmaReal;};

    void makeSigma(        double  & sigma) {sigma = sigmaReal;};

    void makeSigma(complex< float> & sigma) {sigma = complex< float>((float) sigmaReal, (float) sigmaImag);};

    void makeSigma(complex<double> & sigma) {sigma = complex<double>(sigmaReal, sigmaImag);};

    virtual int initSolver(SLV & solver) = 0;

    virtual void dumpAllParameters() = 0;

  // Public members.

  public:

    // Arpack parameters.

    bool symPb; // Symmetric problem.
    a_int nbEV;
    a_int nbCV;
    double tol;
    double sigmaReal, sigmaImag; // Eigen value translation: look for lambda+sigma instead of lambda.
    bool dumpToFile; // Dump resid and v to arpackSolver.*.out files after solve.
    bool restartFromFile; // Restart solve with resid and v values provided in arpackSolver.*.out files.
    string mag; // Magnitude <=> "which" arpack parameter.
    int maxIt;
    bool schur;
    int verbose;

    // Arpack outputs.

    bool stdPb; // Standard or generalized (= not standard).
    StdVecZ val; // Eigen values.
    StdVecEVZ vec; // Eigen vectors.
    int mode;
    int nbIt;
    double imsTime; // Init mode solver time.
    double rciTime; // Reverse communication interface time.

  // Protected members.

  protected:

    a_int nbDim;

  // Private members.

  private:

    RC * resid; // Saved: enable restart from previous solve.
    RC * v;     // Saved: enable restart from previous solve.
};

// Definition of arpackItrSolver class: specialization of arpackSolver using iterative solvers.
// Note: Eigen provides iterative solvers only for sparse matrices.

// RC: Real or Complex.
// FD: Float or Double.
// EM: Eigen Matrix (sparse only, not dense).
// SLV: Solver.
template<typename RC, typename FD, typename EM, typename SLV>
class arpackItrSolver: public arpackSolver<RC, FD, EM, SLV> {
  // Public methods.

  public:

    arpackItrSolver(): arpackSolver<RC, FD, EM, SLV>() {
      slvTol = 1.e-6;
      slvMaxIt = 100;
      slvILUDropTol = 1.;
      slvILUFillFactor = 2;
    };

    void dumpAllParameters() {
      this->dumpParameters();
      if (this->verbose >= 1) {
        cout << endl << "arpackItrSolver:" << endl;
        cout << endl << "slvTol: " << slvTol << endl;
        cout << endl << "slvMaxIt: " << slvMaxIt << endl;
        cout << endl << "slvILUDropTol: " << slvILUDropTol << endl;
        cout << endl << "slvILUFillFactor: " << slvILUFillFactor << endl;
      }
    };

    virtual int initSolver(Eigen::BiCGSTAB<EM> & solver) {
      // Solve with arpack using sparse matrices and iterative solvers.

      solver.setTolerance(slvTol);
      solver.setMaxIterations(slvMaxIt);

      return 0;
    };

    virtual int initSolver(Eigen::ConjugateGradient<EM> & solver) {
      // Solve with arpack using sparse matrices and iterative solvers.

      solver.setTolerance(slvTol);
      solver.setMaxIterations(slvMaxIt);

      return 0;
    };

    virtual int initSolver(Eigen::BiCGSTAB<EM, Eigen::IncompleteLUT<RC>> & solver) {
      // Solve with arpack using sparse matrices and iterative solvers.

      solver.setTolerance(slvTol);
      solver.setMaxIterations(slvMaxIt);
      solver.preconditioner().setDroptol(slvILUDropTol);
      solver.preconditioner().setFillfactor(slvILUFillFactor);

      return 0;
    };

    virtual int initSolver(Eigen::ConjugateGradient<EM, Eigen::Lower|Eigen::Upper, Eigen::IncompleteLUT<RC>> & solver) {
      // Solve with arpack using sparse matrices and iterative solvers.

      solver.setTolerance(slvTol);
      solver.setMaxIterations(slvMaxIt);
      solver.preconditioner().setDroptol(slvILUDropTol);
      solver.preconditioner().setFillfactor(slvILUFillFactor);

      return 0;
    };

  // Public members.

  public:

    // Iterative solvers parameters.

    double slvTol; // Tolerance of the iterative mode solver.
    int    slvMaxIt; // Maximum number of iterations of the iterative mode solver.
    double slvILUDropTol; // Drop tolerance of the ILU preconditioner (if any) of the iterative mode solver.
    int    slvILUFillFactor; // Fill factor of the ILU preconditioner (if any) of the iterative mode solver.
};

// Definition of arpackDrtSolver class: specialization of arpackSolver using direct solvers.
// Note: Eigen provides direct solvers for both sparse and dense matrices.

// RC: Real or Complex.
// FD: Float or Double.
// EM: Eigen Matrix (sparse or dense).
// SLV: Solver.
template<typename RC, typename FD, typename EM, typename SLV>
class arpackDrtSolver: public arpackSolver<RC, FD, EM, SLV> {
  // Public methods.

  public:

    arpackDrtSolver(): arpackSolver<RC, FD, EM, SLV>() {
      slvPvtThd = 1.e-6;
      slvOffset = 0.;
      slvScale = 1.;
    };

    void dumpAllParameters() {
      this->dumpParameters();
      if (this->verbose >= 1) {
        cout << endl << "arpackDrtSolver:" << endl;
        cout << endl << "slvPvtThd: " << slvPvtThd << endl;
        cout << endl << "slvOffset: " << slvOffset << endl;
        cout << endl << "slvScale: " << slvScale << endl;
      }
    };

    virtual int initSolver(Eigen::SparseLU<EM, Eigen::COLAMDOrdering<int>> & solver) {
      // Solve with arpack using sparse matrices and direct solvers.

      solver.setPivotThreshold(slvPvtThd);

      return 0;
    };

    virtual int initSolver(Eigen::SparseQR<EM, Eigen::COLAMDOrdering<int>> & solver) {
      // Solve with arpack using sparse matrices and direct solvers.

      solver.setPivotThreshold(slvPvtThd);

      return 0;
    };

    virtual int initSolver(Eigen::SimplicialLLT<EM, Eigen::Lower, Eigen::COLAMDOrdering<int>> & solver) {
      // Solve with arpack using sparse matrices and direct solvers.

      solver.setShift(slvOffset, slvScale);

      return 0;
    };

    virtual int initSolver(Eigen::SimplicialLDLT<EM, Eigen::Lower, Eigen::COLAMDOrdering<int>> & solver) {
      // Solve with arpack using sparse matrices and direct solvers.

      solver.setShift(slvOffset, slvScale);

      return 0;
    };

    virtual int initSolver(Eigen::LLT<EM> & solver) {
      // Solve with arpack using dense matrices and direct solvers.

      if (this->mode == 1) return 0;

      solver = Eigen::LLT<EM>(this->nbDim);

      return 0;
    };

    virtual int initSolver(Eigen::LDLT<EM> & solver) {
      // Solve with arpack using dense matrices and direct solvers.

      if (this->mode == 1) return 0;

      solver = Eigen::LDLT<EM>(this->nbDim);

      return 0;
    };

    virtual int initSolver(Eigen::FullPivLU<EM> & solver) {
      // Solve with arpack using dense matrices and direct solvers.

      if (this->mode == 1) return 0;

      solver = Eigen::FullPivLU<EM>(this->nbDim, this->nbDim);
      solver.setThreshold(slvPvtThd);

      return 0;
    };

    virtual int initSolver(Eigen::FullPivHouseholderQR<EM> & solver) {
      // Solve with arpack using dense matrices and direct solvers.

      if (this->mode == 1) return 0;

      solver = Eigen::FullPivHouseholderQR<EM>(this->nbDim, this->nbDim);
      solver.setThreshold(slvPvtThd);

      return 0;
    };

    virtual int initSolver(Eigen::PartialPivLU<EM> & solver) {
      // Solve with arpack using dense matrices and direct solvers.

      if (this->mode == 1) return 0;

      solver = Eigen::PartialPivLU<EM>(this->nbDim);

      return 0;
    };

    virtual int initSolver(Eigen::HouseholderQR<EM> & solver) {
      // Solve with arpack using dense matrices and direct solvers.

      if (this->mode == 1) return 0;

      solver = Eigen::HouseholderQR<EM>(this->nbDim, this->nbDim);

      return 0;
    };

  // Public members.

  public:

    // Direct solvers parameters.

    double slvPvtThd; // Pivoting tolerance of the direct mode solver.
    double slvOffset; // Cholesky offset (LLT, LDLT) of the direct mode solver.
    double slvScale; // Cholesky scale (LLT, LDLT) of the direct mode solver.
};

#endif

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=2 ts=2 et smartindent :*/
