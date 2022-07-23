!\BeginDoc
!
!\Name: pdnaupd
!
! Message Passing Layer: MPI
!
!\Description:
!  Reverse communication interface for the Implicitly Restarted Arnoldi
!  iteration. This subroutine computes approximations to a few eigenpairs
!  of a linear operator "OP" with respect to a semi-inner product defined by
!  a symmetric positive semi-definite real matrix B. B may be the identity
!  matrix. NOTE: If the linear operator "OP" is real and symmetric
!  with respect to the real positive semi-definite symmetric matrix B,
!  i.e. B*OP = (OP`)*B, then subroutine dsaupd  should be used instead.
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  pdnaupd  is usually called iteratively to solve one of the
!  following problems:
!
!  Mode 1:  A*x = lambda*x.
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)
!
!  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then
!           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
!           Note: If sigma is real, i.e. imaginary part of sigma is zero;
!                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
!                 amu == 1/(lambda-sigma).
!
!  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then
!           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
!
!  Both mode 3 and 4 give the same enhancement to eigenvalues close to
!  the (complex) shift sigma.  However, as lambda goes to infinity,
!  the operator OP in mode 4 dampens the eigenvalues more strongly than
!  does OP defined in mode 3.
!
!  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
!        should be accomplished either by a direct method
!        using a sparse matrix factorization and solving
!
!           [A - sigma*M]*w = v  or M*w = v,
!
!        or through an iterative method for solving these
!        systems.  If an iterative method is used, the
!        convergence test must be more stringent than
!        the accuracy requirements for the eigenvalue
!        approximations.
!
!\Usage:
!  call pdnaupd
!     ( COMM, IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!\Arguments
!  COMM    MPI  Communicator for the processor grid.  (INPUT)
!
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to pdnaupd .  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          pdnaupd  with the result.  The operand is given in
!          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    In mode 3 and 4, the vector B * X is already
!                    available in WORKD(ipntr(3)).  It does not
!                    need to be recomputed in forming OP * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO =  3: compute the IPARAM(8) real and imaginary parts
!                    of the shifts where INPTR(14) is the pointer
!                    into WORKL for placing the shifts. See Remark
!                    5 below.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
!          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  WHICH   Character*2.  (INPUT)
!          'LM' -> want the NEV eigenvalues of largest magnitude.
!          'SM' -> want the NEV eigenvalues of smallest magnitude.
!          'LR' -> want the NEV eigenvalues of largest real part.
!          'SR' -> want the NEV eigenvalues of smallest real part.
!          'LI' -> want the NEV eigenvalues of largest imaginary part.
!          'SI' -> want the NEV eigenvalues of smallest imaginary part.
!
!  NEV     Integer.  (INPUT)
!          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
!
!  TOL     Double precision  scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
!          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
!          DEFAULT = DLAMCH ('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine DLAMCH ).
!
!  RESID   Double precision  array of length N.  (INPUT/OUTPUT)
!          On INPUT:
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V. NCV must satisfy the two
!          inequalities 2 <= NCV-NEV and NCV <= N.
!          This will indicate how many Arnoldi vectors are generated
!          at each iteration.  After the startup phase in which NEV
!          Arnoldi vectors are generated, the algorithm generates
!          approximately NCV-NEV Arnoldi vectors at each subsequent update
!          iteration. Most of the cost in generating each Arnoldi vector is
!          in the matrix-vector operation OP*x.
!          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz
!          values are kept together. (See remark 4 below)
!
!  V       Double precision  array N by NCV.  (OUTPUT)
!          Contains the final set of Arnoldi basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are provided by the user via
!                      reverse communication.  The real and imaginary
!                      parts of the NCV eigenvalues of the Hessenberg
!                      matrix H are returned in the part of the WORKL
!                      array corresponding to RITZR and RITZI. See remark
!                      5 below.
!          ISHIFT = 1: exact shifts with respect to the current
!                      Hessenberg matrix H.  This is equivalent to
!                      restarting the iteration with a starting vector
!                      that is a linear combination of approximate Schur
!                      vectors associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced.
!
!          IPARAM(3) = MXITER
!          On INPUT:  maximum number of Arnoldi update iterations allowed.
!          On OUTPUT: actual number of Arnoldi update iterations taken.
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" Ritz values.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used.
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4; See under \Description of pdnaupd  for the
!          four modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), pdnaupd  returns NP, the number
!          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
!          5 below.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.
!
!  IPNTR   Integer array of length 14.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD and WORKL
!          arrays for matrices/vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X in WORKD.
!          IPNTR(2): pointer to the current result vector Y in WORKD.
!          IPNTR(3): pointer to the vector B * X in WORKD when used in
!                    the shift-and-invert mode.
!          IPNTR(4): pointer to the next available location in WORKL
!                    that is untouched by the program.
!          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
!                    H in WORKL.
!          IPNTR(6): pointer to the real part of the ritz value array
!                    RITZR in WORKL.
!          IPNTR(7): pointer to the imaginary part of the ritz value array
!                    RITZI in WORKL.
!          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
!                    with the Ritz values located in RITZR and RITZI in WORKL.
!          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
!
!          Note: IPNTR(9:13) is only referenced by dneupd . See Remark 2 below.
!
!          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
!                     original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     pdneupd  if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  WORKD   Double precision  work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration. Upon termination
!          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
!          associated with the converged Ritz values is desired, see remark
!          2 below, subroutine dneupd  uses this output.
!          See Data Distribution Note below.
!
!  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least 3*NCV**2 + 6*NCV.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the
!                Implicitly restarted Arnoldi iteration. One possibility
!                is to increase the size of NCV relative to NEV.
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iteration
!                must be greater than zero.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array is not sufficient.
!          = -8: Error return from LAPACK eigenvalue calculation;
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization.
!
!\Remarks
!  1. The computed Ritz values are approximate eigenvalues of OP. The
!     selection of WHICH should be made with this in mind when
!     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
!     original problem may be obtained with the ARPACK subroutine dneupd .
!
!  2. If a basis for the invariant subspace corresponding to the converged Ritz
!     values is needed, the user must call dneupd  immediately following
!     completion of pdnaupd . This is new starting with release 2 of ARPACK.
!
!  3. If M can be factored into a Cholesky factorization M = LL`
!     then Mode = 2 should not be selected.  Instead one should use
!     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
!     linear systems should be solved with L and L` rather
!     than computing inverses.  After convergence, an approximate
!     eigenvector z of the original problem is recovered by solving
!     L`z = x  where x is a Ritz vector of OP.
!
!  4. At present there is no a-priori analysis to guide the selection
!     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
!     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
!     the same type are to be solved, one should experiment with increasing
!     NCV while keeping NEV fixed for a given test problem.  This will
!     usually decrease the required number of OP*x operations but it
!     also increases the work and storage required to maintain the orthogonal
!     basis vectors.  The optimal "cross-over" with respect to CPU time
!     is problem dependent and must be determined empirically.
!     See Chapter 8 of Reference 2 for further information.
!
!  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
!     NP = IPARAM(8) real and imaginary parts of the shifts in locations
!         real part                  imaginary part
!         -----------------------    --------------
!     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
!     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
!                        .                          .
!                        .                          .
!                        .                          .
!     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
!
!     Only complex conjugate pairs of shifts may be applied and the pairs
!     must be placed in consecutive locations. The real part of the
!     eigenvalues of the current upper Hessenberg matrix are located in
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part
!     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
!     according to the order defined by WHICH. The complex conjugate
!     pairs are kept together and the associated Ritz estimates are located in
!     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------
!
!\Data Distribution Note:
!
!  Fortran-D syntax:
!  ================
!  Double precision  resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!  decompose  d1(n), d2(n,ncv)
!  align      resid(i) with d1(i)
!  align      v(i,j)   with d2(i,j)
!  align      workd(i) with d1(i)     range (1:n)
!  align      workd(i) with d1(i-n)   range (n+1:2*n)
!  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
!  distribute d1(block), d2(block,:)
!  replicated workl(lworkl)
!
!  Cray MPP syntax:
!  ===============
!  Double precision   resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
!  shared     resid(block), v(block,:), workd(block,:)
!  replicated workl(lworkl)
!
!  CM2/CM5 syntax:
!  ==============
!
!-----------------------------------------------------------------------
!
!     include   'ex-nonsym.doc'
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
!     Real Matrices", Linear Algebra and its Applications, vol 88/89,
!     pp 575-595, (1987).
!
!\Routines called:
!     pdnaup2   Parallel ARPACK routine that implements the Implicitly Restarted
!              Arnoldi Iteration.
!     pivout   Parallel ARPACK utility routine that prints integers.
!     arscnd   ARPACK utility routine for timing.
!     pdvout    Parallel ARPACK utility routine that prints vectors.
!     pdlamch10   ScaLAPACK routine that determines machine constants.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              Cray Research, Inc. &
!     Dept. of Computational &     CRPC / Rice University
!     Applied Mathematics          Houston, Texas
!     Rice University
!     Houston, Texas
!
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Code FILE: naupd.F   SID: 2.2
!
!\SCCS Information:
! FILE: naupd.F   SID: 1.9   DATE OF SID: 04/10/01
!
!\Remarks
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pdnaupd
     &   ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &     iparam, ipntr, workd, workl, lworkl, info )
!
      include  'mpif.h'
      include   'pcontext.h'
!
!     %------------------%
!     | MPI Variables    |
!     %------------------%
!
      integer    comm, myid
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &           one, zero
      parameter (one = 1.0 , zero = 0.0 )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    bounds, ierr, ih, iq, ishift, iupd, iw,
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq,
     &           levec, mode, msglvl, mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   pdnaup2 , pdvout , pivout, arscnd, dstatn, pcontext
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &           pdlamch10
      external   pdlamch10
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido .eq. 0) then
!
!        %-------------------------------%
!        | Initialize parallel execution |
!        | context                       |
!        %-------------------------------%
!
      call pcontext
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call dstatn
         call arscnd (t0)
         msglvl = mnaupd
!
!        %----------------%
!        | Error checking |
!        %----------------%
!
         ierr   = 0
         ishift = iparam(1)
!         levec  = iparam(2)
         mxiter = iparam(3)
!         nb     = iparam(4)
         nb     = 1
!
!        %--------------------------------------------%
!        | Revision 2 performs only implicit restart. |
!        %--------------------------------------------%
!
         iupd   = 1
         mode   = iparam(7)
!
         if (n .le. 0) then
             ierr = -1
         else if (nev .le. 0) then
             ierr = -2
         else if (ncv .le. nev+1) then
             ierr = -3
         else if (mxiter .le. 0) then
             ierr = -4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 4) then
                                                ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
                                                ierr = -12
         end if
!
!        %------------%
!        | Error Exit |
!        %------------%
!
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
!
!        %------------------------%
!        | Set default parameters |
!        %------------------------%
!
         if (nb .le. 0)	nb = 1
         if (tol .le. zero)	tol = pdlamch10 (comm, 'EpsMach')
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        | NEV0 is the local variable designating the   |
!        | size of the invariant subspace desired.      |
!        %----------------------------------------------%
!
         np     = ncv - nev
         nev0   = nev
!
!        %-----------------------------%
!        | Zero out internal workspace |
!        %-----------------------------%
!
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
 10      continue
!
!        %-------------------------------------------------------------%
!        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
!        | etc... and the remaining workspace.                         |
!        | Also update pointer to be used on output.                   |
!        | Memory is laid out as follows:                              |
!        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
!        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
!        |                                   parts of ritz values      |
!        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
!        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
!        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
!        | The final workspace is needed by subroutine pdneigh  called  |
!        | by pdnaup2 . Subroutine dneigh  calls LAPACK routines for     |
!        | calculating eigenvalues and the last row of the eigenvector |
!        | matrix.                                                     |
!        %-------------------------------------------------------------%
!
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
!
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw
!
      end if
!
!     %-------------------------------------------------------%
!     | Carry out the Implicitly restarted Arnoldi Iteration. |
!     %-------------------------------------------------------%
!
      call pdnaup2
     &   ( comm, ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr),
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw),
     &     ipntr, workd, info )
!
!     %--------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication |
!     | to compute operations involving OP or shifts.    |
!     %--------------------------------------------------%
!
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
!
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
!
!     %------------------------------------%
!     | Exit if there was an informational |
!     | error within pdnaup2 .              |
!     %------------------------------------%
!
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
!
      if (msglvl .gt. 0) then
         call pivout (comm, logfil, 1, [mxiter], ndigit,
     &               '_naupd: Number of update iterations taken')
         call pivout (comm, logfil, 1, [np], ndigit,
     &               '_naupd: Number of wanted "converged" Ritz values')
         call pdvout  (comm, logfil, np, workl(ritzr), ndigit,
     &               '_naupd: Real part of the final Ritz values')
         call pdvout  (comm, logfil, np, workl(ritzi), ndigit,
     &               '_naupd: Imaginary part of the final Ritz values')
         call pdvout  (comm, logfil, np, workl(bounds), ndigit,
     &               '_naupd: Associated Ritz estimates')
      end if
!
      call arscnd (t1)
      tnaupd = t1 - t0
!
      if (msglvl .gt. 0) then
         call MPI_COMM_RANK( comm, myid, ierr )
         if ( myid .eq. 0 ) then
!
!        %--------------------------------------------------------%
!        | Version Number & Version Date are defined in version.h |
!        %--------------------------------------------------------%
!
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tnaupd, tnaup2, tnaitr, titref,
     &                  tgetv0, tneigh, tngets, tnapps, tnconv, trvec
 1000    format (//,
     &      5x, '=============================================',/
     &      5x, '= Nonsymmetric implicit Arnoldi update code =',/
     &      5x, '= Version Number: ', ' 2.1' , 21x, ' =',/
     &      5x, '= Version Date:   ', ' 3/19/97' , 16x,   ' =',/
     &      5x, '=============================================',/
     &      5x, '= Summary of timing statistics              =',/
     &      5x, '=============================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in p_naup2 routine              = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6,/
     &      5x, 'Total time in computing final Ritz vectors = ', f12.6/)
         end if
      end if
!
 9000 continue
!
      return
!
!     %----------------%
!     | End of pdnaupd  |
!     %----------------%
!
      end
