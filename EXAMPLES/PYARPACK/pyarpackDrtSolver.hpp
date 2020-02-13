#ifndef __PYARPACKDRTSOLVER_HPP__
#define __PYARPACKDRTSOLVER_HPP__

#include <string>

#include <arpackSolver.hpp>
#include <pyarpackServices.hpp>
#include "debug_c.hpp"
#include "stat_c.hpp"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace bp = boost::python;
namespace bn = boost::python::numpy;

template<typename RC, typename FD, typename EM, typename SLV>
class pyarpackSparseDrtSolver: public arpackDrtSolver<RC, FD, EM, SLV> {
  // Public methods.

  public:

    pyarpackSparseDrtSolver(): arpackDrtSolver<RC, FD, EM, SLV>() {
      debug  = 0;

      nopx = 0, nbx = 0, nrorth = 0, nitref = 0, nrstrt = 0;
      tsaupd = 0., tsaup2 = 0., tsaitr = 0., tseigt = 0., tsgets = 0., tsapps = 0., tsconv = 0.;
      tnaupd = 0., tnaup2 = 0., tnaitr = 0., tneigt = 0., tngets = 0., tnapps = 0., tnconv = 0.;
      tcaupd = 0., tcaup2 = 0., tcaitr = 0., tceigt = 0., tcgets = 0., tcapps = 0., tcconv = 0.;
      tmvopx = 0., tmvbx = 0., tgetv0 = 0., titref = 0., trvec = 0.;
    };

    int solve(bp::tuple & A, bp::tuple B = bp::tuple()) {
      ARPACKSOLVERDEBUGSTAT();
      EM M;
      int rc = pyarpackServices<RC, EM>::buildSparseMatrice(A, M, debug, "A");
      if (rc != 0) {pyarpackThrowError("build matrice from A KO"); return rc;}
      bool stdPb = (bp::len(B) > 0) ? false : true;
      EM N;
      if (!stdPb) {
        rc = pyarpackServices<RC, EM>::buildSparseMatrice(B, N, debug, "B");
        if (rc != 0) {pyarpackThrowError("build matrice from B KO"); return rc;}
      }
      return arpackDrtSolver<RC, FD, EM, SLV>::solve(M, (stdPb ? NULL : &N));
    };

    int checkEigVec(bp::tuple const & A, bp::tuple const B = bp::tuple(), double const diffTol = 1.e-3) {
      ARPACKSOLVERDEBUGSTAT();
      EM M;
      int rc = pyarpackServices<RC, EM>::buildSparseMatrice(A, M, debug, "A");
      if (rc != 0) {pyarpackThrowError("build matrice from A KO"); return rc;}
      bool stdPb = (bp::len(B) > 0) ? false : true;
      EM N;
      if (!stdPb) {
        rc = pyarpackServices<RC, EM>::buildSparseMatrice(B, N, debug, "B");
        if (rc != 0) {pyarpackThrowError("build matrice from B KO"); return rc;}
      }
      return arpackDrtSolver<RC, FD, EM, SLV>::checkEigVec(M, (stdPb ? NULL : &N), &diffTol);
    };

  // Public members.

  public:

    a_int debug;

    a_int nopx, nbx, nrorth, nitref, nrstrt;
    float tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv;
    float tnaupd, tnaup2, tnaitr, tneigt, tngets, tnapps, tnconv;
    float tcaupd, tcaup2, tcaitr, tceigt, tcgets, tcapps, tcconv;
    float tmvopx, tmvbx, tgetv0, titref, trvec;
};

template<typename RC, typename FD, typename EM, typename SLV>
class pyarpackDenseDrtSolver: public arpackDrtSolver<RC, FD, EM, SLV> {
  // Public methods.

  public:

    pyarpackDenseDrtSolver(): arpackDrtSolver<RC, FD, EM, SLV>() {
      debug  = 0;

      nopx = 0, nbx = 0, nrorth = 0, nitref = 0, nrstrt = 0;
      tsaupd = 0., tsaup2 = 0., tsaitr = 0., tseigt = 0., tsgets = 0., tsapps = 0., tsconv = 0.;
      tnaupd = 0., tnaup2 = 0., tnaitr = 0., tneigt = 0., tngets = 0., tnapps = 0., tnconv = 0.;
      tcaupd = 0., tcaup2 = 0., tcaitr = 0., tceigt = 0., tcgets = 0., tcapps = 0., tcconv = 0.;
      tmvopx = 0., tmvbx = 0., tgetv0 = 0., titref = 0., trvec = 0.;
    };

    int solve(bp::tuple & A, bp::tuple B = bp::tuple()) {
      ARPACKSOLVERDEBUGSTAT();
      EM M;
      int rc = pyarpackServices<RC, EM>::buildDenseMatrice(A, M, debug, "A");
      if (rc != 0) {pyarpackThrowError("build matrice from A KO"); return rc;}
      bool stdPb = (bp::len(B) > 0) ? false : true;
      EM N;
      if (!stdPb) {
        rc = pyarpackServices<RC, EM>::buildDenseMatrice(B, N, debug, "B");
        if (rc != 0) {pyarpackThrowError("build matrice from B KO"); return rc;}
      }
      return arpackDrtSolver<RC, FD, EM, SLV>::solve(M, (stdPb ? NULL : &N));
    };

    int checkEigVec(bp::tuple const & A, bp::tuple const B = bp::tuple(), double const diffTol = 1.e-3) {
      ARPACKSOLVERDEBUGSTAT();
      EM M;
      int rc = pyarpackServices<RC, EM>::buildDenseMatrice(A, M, debug, "A");
      if (rc != 0) {pyarpackThrowError("build matrice from A KO"); return rc;}
      bool stdPb = (bp::len(B) > 0) ? false : true;
      EM N;
      if (!stdPb) {
        rc = pyarpackServices<RC, EM>::buildDenseMatrice(B, N, debug, "B");
        if (rc != 0) {pyarpackThrowError("build matrice from B KO"); return rc;}
      }
      return arpackDrtSolver<RC, FD, EM, SLV>::checkEigVec(M, (stdPb ? NULL : &N), &diffTol);
    };

  // Public members.

  public:

    a_int debug;

    a_int nopx, nbx, nrorth, nitref, nrstrt;
    float tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv;
    float tnaupd, tnaup2, tnaitr, tneigt, tngets, tnapps, tnconv;
    float tcaupd, tcaup2, tcaitr, tceigt, tcgets, tcapps, tcconv;
    float tmvopx, tmvbx, tgetv0, titref, trvec;
};

#endif

// Local Variables:
// mode: c++
// c-file-style:"stroustrup"
// show-trailing-whitespace: t
// End:
/* vim: set sw=2 ts=2 et smartindent :*/
