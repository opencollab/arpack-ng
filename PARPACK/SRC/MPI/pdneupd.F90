!\BeginDoc
!
!\Name: pdneupd
!
! Message Passing Layer: MPI
!
!\Description:
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  basis is always computed.  There is an additional storage cost of n*nev
!  if both are requested (in this case a separate array Z must be supplied).
!
!  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
!  are derived from approximate eigenvalues and eigenvectors of
!  of the linear operator OP prescribed by the MODE selection in the
!  call to PDNAUPD .  PDNAUPD  must be called before this routine is called.
!  These approximate eigenvalues and vectors are commonly called Ritz
!  values and Ritz vectors respectively.  They are referred to as such
!  in the comments that follow.  The computed orthonormal basis for the
!  invariant subspace corresponding to these Ritz values is referred to as a
!  Schur basis.
!
!  See documentation in the header of the subroutine PDNAUPD  for
!  definition of OP as well as other terms and the relation of computed
!  Ritz values and Ritz vectors of OP with respect to the given problem
!  A*z = lambda*B*z.  For a brief description, see definitions of
!  IPARAM(7), MODE and WHICH in the documentation of PDNAUPD .
!
!\Usage:
!  call pdneupd
!     ( COMM, RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI,
!       WORKEV, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!\Arguments
!  COMM    MPI  Communicator for the processor grid.  (INPUT)
!
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether a basis for the invariant subspace corresponding
!          to the converged Ritz value approximations for the eigenproblem
!          A*z = lambda*B*z is computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
!                                See Remarks below.
!
!  HOWMNY  Character*1  (INPUT)
!          Specifies the form of the basis for the invariant subspace
!          corresponding to the converged Ritz values that is to be computed.
!
!          = 'A': Compute NEV Ritz vectors;
!          = 'P': Compute NEV Schur vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
!
!  DR      Double precision  array of dimension NEV+1.  (OUTPUT)
!          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains
!          the real part of the Ritz  approximations to the eigenvalues of
!          A*z = lambda*B*z.
!          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
!          DR contains the real part of the Ritz values of OP computed by
!          PDNAUPD . A further computation must be performed by the user
!          to transform the Ritz values computed for OP by PDNAUPD  to those
!          of the original system A*z = lambda*B*z. See remark 3 below.
!
!  DI      Double precision  array of dimension NEV+1.  (OUTPUT)
!          On exit, DI contains the imaginary part of the Ritz value
!          approximations to the eigenvalues of A*z = lambda*B*z associated
!          with DR.
!
!          NOTE: When Ritz values are complex, they will come in complex
!                conjugate pairs.  If eigenvectors are requested, the
!                corresponding Ritz vectors will also come in conjugate
!                pairs and the real and imaginary parts of these are
!                represented in two consecutive columns of the array Z
!                (see below).
!
!  Z       Double precision  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
!          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
!          Z represent approximate eigenvectors (Ritz vectors) corresponding
!          to the NCONV=IPARAM(5) Ritz values for eigensystem
!          A*z = lambda*B*z.
!
!          The complex Ritz vector associated with the Ritz value
!          with positive imaginary part is stored in two consecutive
!          columns.  The first column holds the real part of the Ritz
!          vector and the second column holds the imaginary part.  The
!          Ritz vector associated with the Ritz value with negative
!          imaginary part is simply the complex conjugate of the Ritz vector
!          associated with the positive imaginary part.
!
!          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
!
!          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
!          the array Z may be set equal to first NEV+1 columns of the Arnoldi
!          basis array V computed by PDNAUPD .  In this case the Arnoldi basis
!          will be destroyed and overwritten with the eigenvector basis.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
!
!  SIGMAR  Double precision   (INPUT)
!          If IPARAM(7) = 3 or 4, represents the real part of the shift.
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  SIGMAI  Double precision   (INPUT)
!          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift.
!          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
!
!  WORKEV  Double precision  work array of dimension 3*NCV.  (WORKSPACE)
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to PDNAUPD  that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to PDNEUPD  following the last call
!         to PDNAUPD .  These arguments MUST NOT BE MODIFIED between
!         the the last call to PDNAUPD  and the call to PDNEUPD .
!
!  Three of these parameters (V, WORKL, INFO) are also output parameters:
!
!  V       Double precision  N by NCV array.  (INPUT/OUTPUT)
!
!          Upon INPUT: the NCV columns of V contain the Arnoldi basis
!                      vectors for OP as constructed by PDNAUPD  .
!
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
!                       contain approximate Schur vectors that span the
!                       desired invariant subspace.  See Remark 2 below.
!
!          NOTE: If the array Z has been set equal to first NEV+1 columns
!          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
!          Arnoldi basis held by V has been overwritten by the desired
!          Ritz vectors.  If a separate array Z has been passed then
!          the first NCONV=IPARAM(5) columns of V will contain approximate
!          Schur vectors that span the desired invariant subspace.
!
!  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
!          PDNAUPD .  They are not changed by PDNEUPD .
!          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
!          real and imaginary part of the untransformed Ritz values,
!          the upper quasi-triangular matrix for H, and the
!          associated matrix representation of the invariant subspace for H.
!
!          Note: IPNTR(9:13) contains the pointers into WORKL for addresses
!          of the above information computed by PDNEUPD .
!          -------------------------------------------------------------
!          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
!                     original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     PDNEUPD  if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!
!          =  0: Normal exit.
!
!          =  1: The Schur form computed by LAPACK routine dlahqr
!                could not be reordered by LAPACK routine dtrsen .
!                Re-enter subroutine pdneupd  with IPARAM(5)=NCV and
!                increase the size of the arrays DR and DI to have
!                dimension at least dimension NCV and allocate at least NCV
!                columns for Z. NOTE: Not necessary if Z and V share
!                the same space. Please notify the authors if this error
!                occurs.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from calculation of a real Schur form.
!                Informational error from LAPACK routine dlahqr .
!          = -9: Error return from calculation of eigenvectors.
!                Informational error from LAPACK routine dtrevc .
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: HOWMNY = 'S' not yet implemented
!          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
!          = -14: PDNAUPD  did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: PDNEUPD  got a different count of the number of converged
!                 Ritz values than PDNAUPD  got.  This indicates the user
!                 probably made an error in passing data from PDNAUPD  to
!                 PDNEUPD  or that the data was modified before entering
!                 PDNEUPD .
!
!\BeginLib
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
!     pivout  Parallel ARPACK utility routine that prints integers.
!     pdmout   Parallel ARPACK utility routine that prints matrices
!     pdvout   Parallel ARPACK utility routine that prints vectors.
!     dgeqr2   LAPACK routine that computes the QR factorization of
!             a matrix.
!     dlacpy   LAPACK matrix copy routine.
!     dlahqr   LAPACK routine to compute the real Schur form of an
!             upper Hessenberg matrix.
!     pdlamch10  ScaLAPACK routine that determines machine constants.
!     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dlaset   LAPACK matrix initialization routine.
!     dorm2r   LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     dtrevc   LAPACK routine to compute the eigenvectors of a matrix
!             in upper quasi-triangular form.
!     dtrsen   LAPACK routine that re-orders the Schur form.
!     dtrmm    Level 3 BLAS matrix times an upper triangular matrix.
!     dger     Level 2 BLAS rank one update to a matrix.
!     dnrm2    Level 1 BLAS that computes the norm of a vector.
!     dscal    Level 1 BLAS that scales a vector.
!     dcopy    Level 1 BLAS that copies one vector to another .
!
!\Remarks
!
!  1. Currently only HOWMNY = 'A' and 'P' are implemented.
!
!     Let X` denote the transpose of X.
!
!  2. Schur vectors are an orthogonal representation for the basis of
!     Ritz vectors. Thus, their numerical properties are often superior.
!     If RVEC = .TRUE. then the relationship
!             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
!     V(:,1:IPARAM(5))` * V(:,1:IPARAM(5)) = I are approximately satisfied.
!     Here T is the leading submatrix of order IPARAM(5) of the real
!     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
!     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!     each 2-by-2 diagonal block has its diagonal elements equal and its
!     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
!     diagonal block is a complex conjugate pair of Ritz values. The real
!     Ritz values are stored on the diagonal of T.
!
!  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
!     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
!     values computed by PDNAUPD  for OP to those of A*z = lambda*B*z.
!     Set RVEC = .true. and HOWMNY = 'A', and
!     compute
!           Z(:,I)` * A * Z(:,I) if DI(I) = 0.
!     If DI(I) is not equal to zero and DI(I+1) = - D(I),
!     then the desired real and imaginary parts of the Ritz value are
!           Z(:,I)` * A * Z(:,I) +  Z(:,I+1)` * A * Z(:,I+1),
!           Z(:,I)` * A * Z(:,I+1) -  Z(:,I+1)` * A * Z(:,I), respectively.
!     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
!     compute V(:,1:IPARAM(5))` * A * V(:,1:IPARAM(5)) and then an upper
!     quasi-triangular matrix of order IPARAM(5) is computed. See remark
!     2 above.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Code FILE: neupd.F   SID: 2.3
!
!\SCCS Information:
! FILE: neupd.F   SID: 1.8   DATE OF SID: 04/10/01
!
!\EndLib
!
!-----------------------------------------------------------------------
      subroutine pdneupd
     &         (comm , rvec , howmny, select, dr    , di  ,
     &          z    , ldz  , sigmar, sigmai, workev, bmat,
     &          n    , which, nev   , tol   , resid ,
     &          ncv  , v    , ldv   , iparam, ipntr ,
     &          workd, workl, lworkl, info  )
!
!     %--------------------%
!     | MPI  Communicator |
!     %--------------------%
!
      integer   comm
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
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision
     &           sigmar, sigmai, tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           dr(nev+1)    , di(nev+1)    , resid(n)  ,
     &           v(ldv,ncv)   , z(ldz,*)     , workd(3*n),
     &           workl(lworkl), workev(3*ncv)
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
      character  type*6
      integer    bounds, ierr, ih, ihbds,
     &           iheigr, iheigi, iconj , nconv   ,
     &           invsub, iuptri, iwev  , iwork(1),
     &           j     , k     , ldh   , ldq     ,
     &           mode  , msglvl, outncv, ritzr   ,
     &           ritzi , wri   , wrr   , irr     ,
     &           iri   , ibd   , ishift, numcnv  ,
     &           np    , jj
      logical    reord
      Double precision
     &           conds  , rnorm, sep  , temp,
     &           vl(1,1), temp1, eps23
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   dcopy  , dger   , dgeqr2 , dlacpy ,
     &           dlahqr , dlaset , pdmout , dorm2r ,
     &           dtrevc , dtrmm  , dtrsen , dscal  ,
     &           pdvout , pivout
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &           dlapy2 , dnrm2 , pdlamch10
      external   dlapy2 , dnrm2 , pdlamch10
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, min, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
      msglvl = mneupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
      eps23 = pdlamch10 (comm, 'Epsilon-Machine')
      eps23 = eps23**(2.0  / 3.0 )
!
!     %--------------%
!     | Quick return |
!     %--------------%
!
      ierr = 0
!
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
!
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 .and. sigmai .eq. zero) then
         type = 'SHIFTI'
      else if (mode .eq. 3 ) then
         type = 'REALPT'
      else if (mode .eq. 4 ) then
         type = 'IMAGPT'
      else
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
!
!     %--------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
!     | etc... and the remaining workspace.                    |
!     | Also update pointer to be used on output.              |
!     | Memory is laid out as follows:                         |
!     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
!     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
!     |                                   parts of ritz values |
!     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
!     %--------------------------------------------------------%
!
!     %-----------------------------------------------------------%
!     | The following is used and set by PDNEUPD .                 |
!     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
!     |                             real part of the Ritz values. |
!     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
!     |                        imaginary part of the Ritz values. |
!     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
!     |                           error bounds of the Ritz values |
!     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
!     |                             quasi-triangular matrix for H |
!     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
!     |       associated matrix representation of the invariant   |
!     |       subspace for H.                                     |
!     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
!     %-----------------------------------------------------------%
!
      ih     = ipntr(5)
      ritzr  = ipntr(6)
      ritzi  = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheigr = bounds + ldh
      iheigi = iheigr + ldh
      ihbds  = iheigi + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheigr
      ipntr(10) = iheigi
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wrr = 1
      wri = ncv + 1
      iwev = wri + ncv
!
!     %-----------------------------------------%
!     | irr points to the REAL part of the Ritz |
!     |     values computed by _neigh before    |
!     |     exiting _naup2.                     |
!     | iri points to the IMAGINARY part of the |
!     |     Ritz values computed by _neigh      |
!     |     before exiting _naup2.              |
!     | ibd points to the Ritz estimates        |
!     |     computed by _neigh before exiting   |
!     |     _naup2.                             |
!     %-----------------------------------------%
!
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
!
!     %------------------------------------%
!     | RNORM is B-norm of the RESID(1:N). |
!     %------------------------------------%
!
      rnorm = workl(ih+2)
      workl(ih+2) = zero
!
      if (msglvl .gt. 2) then
         call pdvout  (comm, logfil, ncv, workl(irr), ndigit,
     &   '_neupd: Real part of Ritz values passed in from _NAUPD.')
         call pdvout  (comm, logfil, ncv, workl(iri), ndigit,
     &   '_neupd: Imag part of Ritz values passed in from _NAUPD.')
         call pdvout  (comm, logfil, ncv, workl(ibd), ndigit,
     &   '_neupd: Ritz estimates passed in from _NAUPD.')
      end if
!
      if (rvec) then
!
         reord = .false.
!
!        %---------------------------------------------------%
!        | Use the temporary bounds array to store indices   |
!        | These will be used to mark the select array later |
!        %---------------------------------------------------%
!
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
!
!        %-------------------------------------%
!        | Select the wanted Ritz values.      |
!        | Sort the Ritz values so that the    |
!        | wanted ones appear at the tailing   |
!        | NEV positions of workl(irr) and     |
!        | workl(iri).  Move the corresponding |
!        | error estimates in workl(bound)     |
!        | accordingly.                        |
!        %-------------------------------------%
!
         np     = ncv - nev
         ishift = 0
         call pdngets (comm      , ishift       , which     ,
     &                 nev       , np           , workl(irr),
     &                 workl(iri), workl(bounds),
     &                 workl     , workl(np+1))
!
         if (msglvl .gt. 2) then
            call pdvout  (comm, logfil, ncv, workl(irr), ndigit,
     &      '_neupd: Real part of Ritz values after calling _NGETS.')
            call pdvout  (comm, logfil, ncv, workl(iri), ndigit,
     &      '_neupd: Imag part of Ritz values after calling _NGETS.')
            call pdvout  (comm, logfil, ncv, workl(bounds), ndigit,
     &      '_neupd: Ritz value indices after calling _NGETS.')
         end if
!
!        %-----------------------------------------------------%
!        | Record indices of the converged wanted Ritz values  |
!        | Mark the select array for possible reordering       |
!        %-----------------------------------------------------%
!
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23,
     &                 dlapy2  ( workl(irr+ncv-j), workl(iri+ncv-j) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
!
!        %-----------------------------------------------------------%
!        | Check the count (numcnv) of converged Ritz values with    |
!        | the number (nconv) reported by dnaupd.  If these two      |
!        | are different then there has probably been an error       |
!        | caused by incorrect passing of the dnaupd data.           |
!        %-----------------------------------------------------------%
!
         if (msglvl .gt. 2) then
             call pivout(comm, logfil, 1, [numcnv], ndigit,
     &            '_neupd: Number of specified eigenvalues')
             call pivout(comm, logfil, 1, [nconv], ndigit,
     &            '_neupd: Number of "converged" eigenvalues')
         end if
!
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
!
!        %-----------------------------------------------------------%
!        | Call LAPACK routine dlahqr  to compute the real Schur form |
!        | of the upper Hessenberg matrix returned by PDNAUPD .       |
!        | Make a copy of the upper Hessenberg matrix.               |
!        | Initialize the Schur vector matrix Q to the identity.     |
!        %-----------------------------------------------------------%
!
         call dcopy (ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call dlaset ('All', ncv, ncv, zero, one, workl(invsub), ldq)
         call dlahqr (.true.       , .true.       , ncv, 1            ,
     &               ncv          , workl(iuptri), ldh, workl(iheigr),
     &               workl(iheigi), 1            , ncv, workl(invsub),
     &               ldq          , ierr)
         call dcopy  (ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
!
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
!
         if (msglvl .gt. 1) then
            call pdvout  (comm, logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H')
            call pdvout  (comm, logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imaginary part of the Eigenvalues of H')
            call pdvout  (comm, logfil, ncv, workl(ihbds), ndigit,
     &           '_neupd: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call pdmout  (comm, logfil, ncv, ncv,
     &              workl(iuptri), ldh, ndigit,
     &              '_neupd: The upper quasi-triangular matrix ')
            end if
         end if
!
         if (reord) then
!
!           %-----------------------------------------------------%
!           | Reorder the computed upper quasi-triangular matrix. |
!           %-----------------------------------------------------%
!
            call dtrsen ('None'       , 'V'          , select       ,
     &                   ncv          , workl(iuptri), ldh          ,
     &                   workl(invsub), ldq          , workl(iheigr),
     &                   workl(iheigi), nconv        , conds        ,
     &                   sep          , workl(ihbds) , ncv          ,
     &                   iwork        , 1            , ierr         )
!
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
!
            if (msglvl .gt. 2) then
                call pdvout (comm, logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H--reordered')
                call pdvout (comm, logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imag part of the eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call pdmout (comm, logfil, ncv, ncv,
     &             workl(iuptri), ldq, ndigit,
     &             '_neupd: Quasi-triangular matrix after re-ordering')
                end if
            end if
         end if
!
!        %---------------------------------------%
!        | Copy the last row of the Schur vector |
!        | into workl(ihbds).  This will be used |
!        | to compute the Ritz estimates of      |
!        | converged Ritz values.                |
!        %---------------------------------------%
!
         call dcopy (ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
!
!        %----------------------------------------------------%
!        | Place the computed eigenvalues of H into DR and DI |
!        | if a spectral transformation was not used.         |
!        %----------------------------------------------------%
!
         if (type .eq. 'REGULR') then
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
         end if
!
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(invsub,ldq).                            |
!        %----------------------------------------------------------%
!
         call dgeqr2 (ncv, nconv , workl(invsub),
     &               ldq, workev, workev(ncv+1),
     &               ierr)
!
!        %---------------------------------------------------------%
!        | * Postmultiply V by Q using dorm2r .                     |
!        | * Copy the first NCONV columns of VQ into Z.            |
!        | * Postmultiply Z by R.                                  |
!        | The N by NCONV matrix Z is now a matrix representation  |
!        | of the approximate invariant subspace associated with   |
!        | the Ritz values in workl(iheigr) and workl(iheigi)      |
!        | The first NCONV columns of V are now approximate Schur  |
!        | vectors associated with the real upper quasi-triangular |
!        | matrix of order NCONV in workl(iuptri)                  |
!        %---------------------------------------------------------%
!
         call dorm2r ('Right', 'Notranspose', n            ,
     &                ncv    , nconv        , workl(invsub),
     &                ldq    , workev       , v            ,
     &                ldv    , workd(n+1)   , ierr         )
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
!
         do 20 j=1, nconv
!
!           %---------------------------------------------------%
!           | Perform both a column and row scaling if the      |
!           | diagonal element of workl(invsub,ldq) is negative |
!           | I'm lazy and don't take advantage of the upper    |
!           | quasi-triangular form of workl(iuptri,ldq)        |
!           | Note that since Q is orthogonal, R is a diagonal  |
!           | matrix consisting of plus or minus ones           |
!           %---------------------------------------------------%
!
            if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
               call dscal (nconv, -one, workl(iuptri+j-1), ldq)
               call dscal (nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
!
 20      continue
!
         if (howmny .eq. 'A') then
!
!           %--------------------------------------------%
!           | Compute the NCONV wanted eigenvectors of T |
!           | located in workl(iuptri,ldq).              |
!           %--------------------------------------------%
!
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
!
            call dtrevc ('Right', 'Select'     , select       ,
     &                   ncv    , workl(iuptri), ldq          ,
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , ierr)
!
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
!
!           %------------------------------------------------%
!           | Scale the returning eigenvectors so that their |
!           | Euclidean norms are all one. LAPACK subroutine |
!           | dtrevc  returns each eigenvector normalized so  |
!           | that the element of largest magnitude has      |
!           | magnitude 1;                                   |
!           %------------------------------------------------%
!
            iconj = 0
            do 40 j=1, nconv
!
               if ( workl(iheigi+j-1) .eq. zero ) then
!
!                 %----------------------%
!                 | real eigenvalue case |
!                 %----------------------%
!
                  temp = dnrm2 ( ncv, workl(invsub+(j-1)*ldq), 1 )
                  call dscal  ( ncv, one / temp,
     &                 workl(invsub+(j-1)*ldq), 1 )
               else
!
!                 %-------------------------------------------%
!                 | Complex conjugate pair case. Note that    |
!                 | since the real and imaginary part of      |
!                 | the eigenvector are stored in consecutive |
!                 | columns, we further normalize by the      |
!                 | square root of two.                       |
!                 %-------------------------------------------%
!
                  if (iconj .eq. 0) then
                     temp = dlapy2 (dnrm2 (ncv,
     &                                   workl(invsub+(j-1)*ldq),
     &                                   1 ),
     &                             dnrm2 (ncv,
     &                                   workl(invsub+j*ldq),
     &                                   1)
     &                             )
                     call dscal (ncv, one/temp,
     &                          workl(invsub+(j-1)*ldq), 1)
                     call dscal (ncv, one/temp,
     &                          workl(invsub+j*ldq), 1)
                     iconj = 1
                  else
                     iconj = 0
                  end if
!
               end if
!
 40         continue
!
            call dgemv ('T'         , ncv          , nconv,
     &                 one         , workl(invsub), ldq  ,
     &                 workl(ihbds), 1            , zero ,
     &                 workev      , 1)
!
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) .ne. zero) then
!
!                 %-------------------------------------------%
!                 | Complex conjugate pair case. Note that    |
!                 | since the real and imaginary part of      |
!                 | the eigenvector are stored in consecutive |
!                 %-------------------------------------------%
!
                  if (iconj .eq. 0) then
                     workev(j) = dlapy2 (workev(j), workev(j+1))
                     workev(j+1) = workev(j)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               end if
 45         continue
!
            if (msglvl .gt. 2) then
               call pdvout (comm, logfil, ncv, workl(ihbds), ndigit,
     &              '_neupd: Last row of the eigenvector matrix for T')
               if (msglvl .gt. 3) then
                  call pdmout (comm, logfil, ncv, ncv,
     &               workl(invsub), ldq, ndigit,
     &               '_neupd: The eigenvector matrix for T')
               end if
            end if
!
!
!           %---------------------------------------%
!           | Copy Ritz estimates into workl(ihbds) |
!           %---------------------------------------%
!
            call dcopy (nconv, workev, 1, workl(ihbds), 1)
!
!           %---------------------------------------------------------%
!           | Compute the QR factorization of the eigenvector matrix  |
!           | associated with leading portion of T in the first NCONV |
!           | columns of workl(invsub,ldq).                           |
!           %---------------------------------------------------------%
!
            call dgeqr2 (ncv, nconv , workl(invsub),
     &                   ldq, workev, workev(ncv+1),
     &                   ierr)
!
!           %----------------------------------------------%
!           | * Postmultiply Z by Q.                       |
!           | * Postmultiply Z by R.                       |
!           | The N by NCONV matrix Z is now contains the  |
!           | Ritz vectors associated with the Ritz values |
!           | in workl(iheigr) and workl(iheigi).          |
!           %----------------------------------------------%
!
            call dorm2r ('Right', 'Notranspose', n            ,
     &                   ncv    , nconv        , workl(invsub),
     &                   ldq    , workev       , z            ,
     &                   ldz    , workd(n+1)   , ierr)
!
            call dtrmm ('Right'   , 'Upper'      , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
!
         end if
!
      else
!
!        %------------------------------------------------------%
!        | An approximate invariant subspace is not needed.     |
!        | Place the Ritz values computed PDNAUPD  into DR and DI |
!        %------------------------------------------------------%
!
         call dcopy (nconv, workl(ritzr), 1, dr, 1)
         call dcopy (nconv, workl(ritzi), 1, di, 1)
         call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
         call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
         call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
      end if
!
!     %------------------------------------------------%
!     | Transform the Ritz values and possibly vectors |
!     | and corresponding error bounds of OP to those  |
!     | of A*x = lambda*B*x.                           |
!     %------------------------------------------------%
!
      if (type .eq. 'REGULR') then
!
         if (rvec)
     &      call dscal (ncv, rnorm, workl(ihbds), 1)
!
      else
!
!        %---------------------------------------%
!        |   A spectral transformation was used. |
!        | * Determine the Ritz estimates of the |
!        |   Ritz values in the original system. |
!        %---------------------------------------%
!
         if (type .eq. 'SHIFTI') then
!
            if (rvec)
     &         call dscal (ncv, rnorm, workl(ihbds), 1)
            do 50 k=1, ncv
               temp = dlapy2 (workl(iheigr+k-1),
     &                       workl(iheigi+k-1) )
               workl(ihbds+k-1) = abs( workl(ihbds+k-1) )
     &                          / temp / temp
 50         continue
!
         else if (type .eq. 'REALPT') then
!
            do 60 k=1, ncv
 60         continue
!
         else if (type .eq. 'IMAGPT') then
!
            do 70 k=1, ncv
 70         continue
!
         end if
!
!        %-----------------------------------------------------------%
!        | *  Transform the Ritz values back to the original system. |
!        |    For TYPE = 'SHIFTI' the transformation is              |
!        |             lambda = 1/theta + sigma                      |
!        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
!        |    Rayleigh quotients or a projection. See remark 3 above.|
!        | NOTES:                                                    |
!        | *The Ritz vectors are not affected by the transformation. |
!        %-----------------------------------------------------------%
!
         if (type .eq. 'SHIFTI') then
!
            do 80 k=1, ncv
               temp = dlapy2 (workl(iheigr+k-1),
     &                       workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1) / temp / temp
     &                           + sigmar
               workl(iheigi+k-1) = -workl(iheigi+k-1) / temp / temp
     &                           + sigmai
 80         continue
!
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
!
         else if (type .eq. 'REALPT' .or. type .eq. 'IMAGPT') then
!
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
!
         end if
!
      if (type .eq. 'SHIFTI' .and. msglvl .gt. 1) then
         call pdvout  (comm, logfil, nconv, dr, ndigit,
     &   '_neupd: Untransformed real part of the Ritz valuess.')
         call pdvout  (comm, logfil, nconv, di, ndigit,
     &   '_neupd: Untransformed imag part of the Ritz valuess.')
         call pdvout  (comm, logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Ritz estimates of untransformed Ritz values.')
      else if (type .eq. 'REGULR' .and. msglvl .gt. 1) then
         call pdvout  (comm, logfil, nconv, dr, ndigit,
     &   '_neupd: Real parts of converged Ritz values.')
         call pdvout  (comm, logfil, nconv, di, ndigit,
     &   '_neupd: Imag parts of converged Ritz values.')
         call pdvout  (comm, logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Associated Ritz estimates.')
      end if
!
      end if
!
!     %-------------------------------------------------%
!     | Eigenvector Purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 2.                                   |
!     %-------------------------------------------------%
!
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
!
!        %------------------------------------------------%
!        | Purify the computed Ritz vectors by adding a   |
!        | little bit of the residual vector:             |
!        |                      T                         |
!        |          resid(:)*( e    s ) / theta           |
!        |                      NCV                       |
!        | where H s = s theta. Remember that when theta  |
!        | has nonzero imaginary part, the corresponding  |
!        | Ritz vector is stored across two columns of Z. |
!        %------------------------------------------------%
!
         iconj = 0
         do 110 j=1, nconv
            if (workl(iheigi+j-1) .eq. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1)
     &                   /  workl(iheigr+j-1)
            else if (iconj .eq. 0) then
               temp = dlapy2 ( workl(iheigr+j-1), workl(iheigi+j-1) )
               workev(j) = ( workl(invsub+(j-1)*ldq+ncv-1) *
     &                       workl(iheigr+j-1) +
     &                       workl(invsub+j*ldq+ncv-1) *
     &                       workl(iheigi+j-1) ) / temp / temp
               workev(j+1) = ( workl(invsub+j*ldq+ncv-1) *
     &                         workl(iheigr+j-1) -
     &                         workl(invsub+(j-1)*ldq+ncv-1) *
     &                         workl(iheigi+j-1) ) / temp / temp
               iconj = 1
            else
               iconj = 0
            end if
 110     continue
!
!        %---------------------------------------%
!        | Perform a rank one update to Z and    |
!        | purify all the Ritz vectors together. |
!        %---------------------------------------%
!
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
!
      end if
!
 9000 continue
!
      return
!
!     %----------------%
!     | End of PDNEUPD  |
!     %----------------%
!
      end
