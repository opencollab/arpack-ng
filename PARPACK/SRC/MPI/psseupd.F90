!\BeginDoc
!
!\Name: psseupd
!
! Message Passing Layer: MPI
!
!\Description:
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) the corresponding approximate eigenvectors,
!
!      (2) an orthonormal (Lanczos) basis for the associated approximate
!          invariant subspace,
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  (Lanczos) basis is always computed.  There is an additional storage cost
!  of n*nev if both are requested (in this case a separate array Z must be
!  supplied).
!
!  These quantities are obtained from the Lanczos factorization computed
!  by PSSAUPD for the linear operator OP prescribed by the MODE selection
!  (see IPARAM(7) in PSSAUPD documentation.)  PSSAUPD must be called before
!  this routine is called. These approximate eigenvalues and vectors are
!  commonly called Ritz values and Ritz vectors respectively.  They are
!  referred to as such in the comments that follow.   The computed orthonormal
!  basis for the invariant subspace corresponding to these Ritz values is
!  referred to as a Lanczos basis.
!
!  See documentation in the header of the subroutine PSSAUPD for a definition
!  of OP as well as other terms and the relation of computed Ritz values
!  and vectors of OP with respect to the given problem  A*z = lambda*B*z.
!
!  The approximate eigenvalues of the original problem are returned in
!  ascending algebraic order.  The user may elect to call this routine
!  once for each desired Ritz vector and store it peripherally if desired.
!  There is also the option of computing a selected set of these vectors
!  with a single call.
!
!\Usage:
!  call psseupd
!     ( COMM, RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
!       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether Ritz vectors corresponding to the Ritz value
!          approximations to the eigenproblem A*z = lambda*B*z are computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute Ritz vectors.
!
!  HOWMNY  Character*1  (INPUT)
!          Specifies how many Ritz vectors are wanted and the form of Z
!          the matrix of Ritz vectors. See remark 1 below.
!          = 'A': compute NEV Ritz vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A' , SELECT is used as workspace.
!
!  D       Real array of dimension NEV.  (OUTPUT)
!          On exit, D contains the Ritz value approximations to the
!          eigenvalues of A*z = lambda*B*z. The values are returned
!          in ascending order. If IPARAM(7) = 3,4,5 then D represents
!          the Ritz values of OP computed by pssaupd transformed to
!          those of the original eigensystem A*z = lambda*B*z. If
!          IPARAM(7) = 1,2 then the Ritz values of OP are the same
!          as the those of A*z = lambda*B*z.
!
!  Z       Real N by NEV array if HOWMNY = 'A'.  (OUTPUT)
!          On exit, Z contains the B-orthonormal Ritz vectors of the
!          eigensystem A*z = lambda*B*z corresponding to the Ritz
!          value approximations.
!          If  RVEC = .FALSE. then Z is not referenced.
!          NOTE: The array Z may be set equal to first NEV columns of the
!          Arnoldi/Lanczos basis array V computed by PSSAUPD.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
!
!  SIGMA   Real  (INPUT)
!          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
!          IPARAM(7) = 1 or 2.
!
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to PSNAUPD that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to PSSEUPD following the last call
!         to PSSAUPD.  These arguments MUST NOT BE MODIFIED between
!         the the last call to PSSAUPD and the call to PSSEUPD.
!
!  Two of these parameters (WORKL, INFO) are also output parameters:
!
!  WORKL   Real work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:4*ncv) contains information obtained in
!          PSSAUPD. They are not changed by PSSEUPD.
!          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
!          untransformed Ritz values, the computed error estimates,
!          and the associated eigenvector matrix of H.
!
!          Note: IPNTR(8:10) contains the pointers into WORKL for addresses
!          of the above information computed by PSSEUPD.
!          -------------------------------------------------------------
!          IPNTR(8): pointer to the NCV RITZ values of the original system.
!          IPNTR(9): pointer to the NCV corresponding error bounds.
!          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
!                     of the tridiagonal matrix T. Only referenced by
!                     PSSEUPD if RVEC = .TRUE. See Remarks.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV must be greater than NEV and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Information error from LAPACK routine ssteqr.
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: NEV and WHICH = 'BE' are incompatible.
!          = -14: PSSAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
!          = -16: HOWMNY = 'S' not yet implemented
!          = -17: SSEUPD got a different count of the number of converged
!                 Ritz values than SSAUPD got.  This indicates the user
!                 probably made an error in passing data from SSAUPD to
!                 SSEUPD or that the data was modified before entering
!                 SSEUPD.
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
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Remarks
!  1. The converged Ritz values are always returned in increasing
!     (algebraic) order.
!
!  2. Currently only HOWMNY = 'A' is implemented. It is included at this
!     stage for the user who wants to incorporate it.
!
!\Routines called:
!     ssesrt  ARPACK routine that sorts an array X, and applies the
!             corresponding permutation to a matrix A.
!     ssortr  ssortr  ARPACK sorting routine.
!     psnorm2 Parallel ARPACK routine that computes the 2-norm of a vector.
!     pivout  Parallel ARPACK utility routine that prints integers.
!     psvout  Parallel ARPACK utility routine that prints vectors.
!     sgeqr2  LAPACK routine that computes the QR factorization of
!             a matrix.
!     slacpy  LAPACK matrix copy routine.
!     pslamch10 ScaLAPACK routine that determines machine constants.
!     sorm2r  LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     ssteqr  LAPACK routine that computes eigenvalues and eigenvectors
!             of a tridiagonal matrix.
!     sger    Level 2 BLAS rank one update to a matrix.
!     scopy   Level 1 BLAS that copies one vector to another .
!     sscal   Level 1 BLAS that scales a vector.
!     sswap   Level 1 BLAS that swaps the contents of two vectors.
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
!     Starting Point: Serial Code FILE: seupd.F   SID: 2.4
!
!\SCCS Information:
! FILE: seupd.F   SID: 1.11   DATE OF SID: 10/25/03
!
!\EndLib
!
!-----------------------------------------------------------------------
      subroutine psseupd
     &    (comm  , rvec  , howmny, select, d    ,
     &     z     , ldz   , sigma , bmat  , n    ,
     &     which , nev   , tol   , resid , ncv  ,
     &     v     , ldv   , iparam, ipntr , workd,
     &     workl , lworkl, info )
!
!     %--------------------%
!     | MPI Communicator |
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
      Real
     &           sigma, tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Real
     &           d(nev), resid(n), v(ldv,ncv), z(ldz, nev),
     &           workd(2*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  type*6
      integer    bounds , ierr   , ih    , ihb   , ihd   ,
     &           iq     , iw     , j     , k     , ldh   ,
     &           ldq    , mode   , msglvl, nconv , next  ,
     &           ritz   , irz    , ibd   , np    , ishift,
     &           leftptr, rghtptr, numcnv, jj
      Real
     &           bnorm2, rnorm, temp, temp1, eps23
      logical    reord
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy , sger  , sgeqr2, slacpy, sorm2r, sscal,
     &           ssesrt, ssteqr, sswap , psvout, pivout, ssortr
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           psnorm2, pslamch10
      external   psnorm2, pslamch10
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
      msglvl = mseupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
!
!     %--------------%
!     | Quick return |
!     %--------------%
!
      if (nconv .eq. 0) go to 9000
      ierr = 0
!
      if (nconv .le. 0)                        ierr = -14
      if (n .le. 0)                            ierr = -1
      if (nev .le. 0)                          ierr = -2
      if (ncv .le. nev)                        ierr = -3
      if (which .ne. 'LM' .and.
     &    which .ne. 'SM' .and.
     &    which .ne. 'LA' .and.
     &    which .ne. 'SA' .and.
     &    which .ne. 'BE')                     ierr = -5
      if (bmat .ne. 'I' .and. bmat .ne. 'G')   ierr = -6
      if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec )
     &                                         ierr = -15
      if (rvec .and. howmny .eq. 'S')           ierr = -16
!
      if (rvec .and. lworkl .lt. ncv**2+8*ncv) ierr = -7
!
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else if (mode .eq. 4 ) then
         type = 'BUCKLE'
      else if (mode .eq. 5 ) then
         type = 'CAYLEY'
      else
                                               ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')     ierr = -11
      if (nev .eq. 1 .and. which .eq. 'BE')    ierr = -12
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
!     %-------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
!     | etc... and the remaining workspace.                   |
!     | Also update pointer to be used on output.             |
!     | Memory is laid out as follows:                        |
!     | workl(1:2*ncv) := generated tridiagonal matrix H      |
!     |       The subdiagonal is stored in workl(2:ncv).      |
!     |       The dead spot is workl(1) but upon exiting      |
!     |       pssaupd stores the B-norm of the last residual  |
!     |       vector in workl(1). We use this !!!             |
!     | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
!     |       The wanted values are in the first NCONV spots. |
!     | workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   |
!     |       The wanted values are in the first NCONV spots. |
!     | NOTE: workl(1:4*ncv) is set by pssaupd and is not     |
!     |       modified by psseupd.                            |
!     %-------------------------------------------------------%
!
!     %-------------------------------------------------------%
!     | The following is used and set by psseupd.             |
!     | workl(4*ncv+1:4*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the diagonal of H. Upon EXIT contains the NCV   |
!     |       Ritz values of the original system. The first   |
!     |       NCONV spots have the wanted values. If MODE =   |
!     |       1 or 2 then will equal workl(2*ncv+1:3*ncv).    |
!     | workl(5*ncv+1:5*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the subdiagonal of H. Upon EXIT contains the    |
!     |       NCV corresponding Ritz estimates of the         |
!     |       original system. The first NCONV spots have the |
!     |       wanted values. If MODE = 1,2 then will equal    |
!     |       workl(3*ncv+1:4*ncv).                           |
!     | workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  |
!     |       the eigenvector matrix for H as returned by     |
!     |       ssteqr. Not referenced if RVEC = .False.        |
!     |       Ordering follows that of workl(4*ncv+1:5*ncv)   |
!     | workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         |
!     |       Workspace. Needed by ssteqr and by psseupd.    |
!     | GRAND total of NCV*(NCV+8) locations.                 |
!     %-------------------------------------------------------%
!
!
      ih     = ipntr(5)
      ritz   = ipntr(6)
      bounds = ipntr(7)
      ldh    = ncv
      ldq    = ncv
      ihd    = bounds + ldh
      ihb    = ihd    + ldh
      iq     = ihb    + ldh
      iw     = iq     + ldh*ncv
      next   = iw     + 2*ncv
      ipntr(4)  = next
      ipntr(8)  = ihd
      ipntr(9)  = ihb
      ipntr(10) = iq
!
!     %----------------------------------------%
!     | irz points to the Ritz values computed |
!     |     by _seigt before exiting _saup2.   |
!     | ibd points to the Ritz estimates       |
!     |     computed by _seigt before exiting  |
!     |     _saup2.                            |
!     %----------------------------------------%
!
      irz = ipntr(11)+ncv
      ibd = irz+ncv
!
!
!     %---------------------------------%
!     | Set machine dependent constant. |
!     %---------------------------------%
!
      eps23 = pslamch10(comm, 'Epsilon-Machine')
      eps23 = eps23**(2.0 / 3.0)
!
!     %---------------------------------------%
!     | RNORM is B-norm of the RESID(1:N).    |
!     | BNORM2 is the 2 norm of B*RESID(1:N). |
!     | Upon exit of pssaupd WORKD(1:N) has   |
!     | B*RESID(1:N).                         |
!     %---------------------------------------%
!
      rnorm = workl(ih)
      if (bmat .eq. 'I') then
         bnorm2 = rnorm
      else if (bmat .eq. 'G') then
         bnorm2 = psnorm2(comm, n, workd, 1)
      end if
!
      if (msglvl .gt. 2) then
         call psvout(comm, logfil, ncv, workl(irz), ndigit,
     &   '_seupd: Ritz values passed in from _SAUPD.')
         call psvout(comm, logfil, ncv, workl(ibd), ndigit,
     &   '_seupd: Ritz estimates passed in from _SAUPD.')
      end if
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
!        | NEV positions of workl(irz).  Move  |
!        | the corresponding error estimates   |
!        | in workl(bound) accordingly.        |
!        %-------------------------------------%
!
         np     = ncv - nev
         ishift = 0
         call pssgets(comm         , ishift, which     ,
     &                nev          , np    , workl(irz),
     &                workl(bounds), workl)
!
         if (msglvl .gt. 2) then
            call psvout(comm, logfil, ncv, workl(irz), ndigit,
     &      '_seupd: Ritz values after calling _SGETS.')
            call psvout(comm, logfil, ncv, workl(bounds), ndigit,
     &      '_seupd: Ritz value indices after calling _SGETS.')
         end if
!
!        %-----------------------------------------------------%
!        | Record indices of the converged wanted Ritz values  |
!        | Mark the select array for possible reordering       |
!        %-----------------------------------------------------%
!
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23, abs(workl(irz+ncv-j)) )
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
!        | the number (nconv) reported by _saupd.  If these two      |
!        | are different then there has probably been an error       |
!        | caused by incorrect passing of the _saupd data.           |
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
            info = -17
            go to 9000
         end if
!
!        %-----------------------------------------------------------%
!        | Call LAPACK routine _steqr to compute the eigenvalues and |
!        | eigenvectors of the final symmetric tridiagonal matrix H. |
!        | Initialize the eigenvector matrix Q to the identity.      |
!        %-----------------------------------------------------------%
!
         call scopy (ncv-1, workl(ih+1)  , 1, workl(ihb), 1)
         call scopy (ncv  , workl(ih+ldh), 1, workl(ihd), 1)
!
         call ssteqr('Identity', ncv      , workl(ihd),
     &               workl(ihb), workl(iq), ldq       ,
     &               workl(iw) , ierr)
!
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
!
         if (msglvl .gt. 1) then
            call scopy (ncv, workl(iq+ncv-1), ldq, workl(iw), 1)
            call psvout (comm, logfil, ncv, workl(ihd), ndigit,
     &          '_seupd: NCV Ritz values of the final H matrix')
            call psvout (comm, logfil, ncv, workl(iw), ndigit,
     &           '_seupd: last row of the eigenvector matrix for H')
         end if
!
         if (reord) then
!
!           %---------------------------------------------%
!           | Reordered the eigenvalues and eigenvectors  |
!           | computed by _steqr so that the "converged"  |
!           | eigenvalues appear in the first NCONV       |
!           | positions of workl(ihd), and the associated |
!           | eigenvectors appear in the first NCONV      |
!           | columns.                                    |
!           %---------------------------------------------%
!
            leftptr = 1
            rghtptr = ncv
!
            if (ncv .eq. 1) go to 30
!
 20         if (select(leftptr)) then
!
!              %-------------------------------------------%
!              | Search, from the left, for the first Ritz |
!              | value that has not converged.             |
!              %-------------------------------------------%
!
               leftptr = leftptr + 1
!
            else if (.not.select(rghtptr)) then
!
!              %----------------------------------------------%
!              | Search, from the right, the first Ritz value |
!              | that has converged.                          |
!              %----------------------------------------------%
!
               rghtptr = rghtptr - 1
!
            else
!
!              %----------------------------------------------%
!              | Swap the Ritz value on the left that has not |
!              | converged with the Ritz value on the right   |
!              | that has converged.  Swap the associated     |
!              | eigenvector of the tridiagonal matrix H as   |
!              | well.                                        |
!              %----------------------------------------------%
!
               temp = workl(ihd+leftptr-1)
               workl(ihd+leftptr-1) = workl(ihd+rghtptr-1)
               workl(ihd+rghtptr-1) = temp
               call scopy(ncv, workl(iq+ncv*(leftptr-1)), 1,
     &                    workl(iw), 1)
               call scopy(ncv, workl(iq+ncv*(rghtptr-1)), 1,
     &                    workl(iq+ncv*(leftptr-1)), 1)
               call scopy(ncv, workl(iw), 1,
     &                    workl(iq+ncv*(rghtptr-1)), 1)
               leftptr = leftptr + 1
               rghtptr = rghtptr - 1
!
            end if
!
            if (leftptr .lt. rghtptr) go to 20
!
 30      end if
!
         if (msglvl .gt. 2) then
             call psvout (comm, logfil, ncv, workl(ihd), ndigit,
     &       '_seupd: The eigenvalues of H--reordered')
         end if
!
!        %----------------------------------------%
!        | Load the converged Ritz values into D. |
!        %----------------------------------------%
!
         call scopy(nconv, workl(ihd), 1, d, 1)
!
      else
!
!        %-----------------------------------------------------%
!        | Ritz vectors not required. Load Ritz values into D. |
!        %-----------------------------------------------------%
!
         call scopy(nconv, workl(ritz), 1, d, 1)
         call scopy(ncv, workl(ritz), 1, workl(ihd), 1)
!
      end if
!
!     %------------------------------------------------------------------%
!     | Transform the Ritz values and possibly vectors and corresponding |
!     | Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values |
!     | (and corresponding data) are returned in ascending order.        |
!     %------------------------------------------------------------------%
!
      if (type .eq. 'REGULR') then
!
!        %---------------------------------------------------------%
!        | Ascending sort of wanted Ritz values, vectors and error |
!        | bounds. Not necessary if only Ritz values are desired.  |
!        %---------------------------------------------------------%
!
         if (rvec) then
            call ssesrt('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call scopy(ncv, workl(bounds), 1, workl(ihb), 1)
         end if
!
      else
!
!        %-------------------------------------------------------------%
!        | *  Make a copy of all the Ritz values.                      |
!        | *  Transform the Ritz values back to the original system.   |
!        |    For TYPE = 'SHIFTI' the transformation is                |
!        |             lambda = 1/theta + sigma                        |
!        |    For TYPE = 'BUCKLE' the transformation is                |
!        |             lambda = sigma * theta / ( theta - 1 )          |
!        |    For TYPE = 'CAYLEY' the transformation is                |
!        |             lambda = sigma * (theta + 1) / (theta - 1 )     |
!        |    where the theta are the Ritz values returned by pssaupd. |
!        | NOTES:                                                      |
!        | *The Ritz vectors are not affected by the transformation.   |
!        |  They are only reordered.                                   |
!        %-------------------------------------------------------------%
!
         call scopy (ncv, workl(ihd), 1, workl(iw), 1)
         if (type .eq. 'SHIFTI') then
            do 40 k=1, ncv
               workl(ihd+k-1) = one / workl(ihd+k-1) + sigma
  40        continue
         else if (type .eq. 'BUCKLE') then
            do 50 k=1, ncv
               workl(ihd+k-1) = sigma * workl(ihd+k-1) /
     &                          (workl(ihd+k-1) - one)
  50        continue
         else if (type .eq. 'CAYLEY') then
            do 60 k=1, ncv
               workl(ihd+k-1) = sigma * (workl(ihd+k-1) + one) /
     &                          (workl(ihd+k-1) - one)
  60        continue
         end if
!
!        %-------------------------------------------------------------%
!        | *  Store the wanted NCONV lambda values into D.             |
!        | *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   |
!        |    into ascending order and apply sort to the NCONV theta   |
!        |    values in the transformed system. We will need this to   |
!        |    compute Ritz estimates in the original system.           |
!        | *  Finally sort the lambda`s into ascending order and apply |
!        |    to Ritz vectors if wanted. Else just sort lambda`s into  |
!        |    ascending order.                                         |
!        | NOTES:                                                      |
!        | *workl(iw:iw+ncv-1) contain the theta ordered so that they  |
!        |  match the ordering of the lambda. We`ll use them again for |
!        |  Ritz vector purification.                                  |
!        %-------------------------------------------------------------%
!
         call scopy (nconv, workl(ihd), 1, d, 1)
         call ssortr('LA', .true., nconv, workl(ihd), workl(iw))
         if (rvec) then
            call ssesrt('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call scopy(ncv, workl(bounds), 1, workl(ihb), 1)
            call sscal(ncv, bnorm2/rnorm, workl(ihb), 1)
            call ssortr('LA', .true., nconv, d, workl(ihb))
         end if
!
      end if
!
!     %------------------------------------------------%
!     | Compute the Ritz vectors. Transform the wanted |
!     | eigenvectors of the symmetric tridiagonal H by |
!     | the Lanczos basis matrix V.                    |
!     %------------------------------------------------%
!
      if (rvec .and. howmny .eq. 'A') then
!
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(iq,ldq).                                |
!        %----------------------------------------------------------%
!
         call sgeqr2(ncv, nconv        , workl(iq) ,
     &               ldq, workl(iw+ncv), workl(ihb),
     &               ierr)
!
!        %--------------------------------------------------------%
!        | * Postmultiply V by Q.                                 |
!        | * Copy the first NCONV columns of VQ into Z.           |
!        | The N by NCONV matrix Z is now a matrix representation |
!        | of the approximate invariant subspace associated with  |
!        | the Ritz values in workl(ihd).                         |
!        %--------------------------------------------------------%
!
         call sorm2r('Right'      , 'Notranspose', n        ,
     &                ncv          , nconv        , workl(iq),
     &                ldq          , workl(iw+ncv), v        ,
     &                ldv          , workd(n+1)   , ierr     )
         call slacpy('All', n, nconv, v, ldv, z, ldz)
!
!        %-----------------------------------------------------%
!        | In order to compute the Ritz estimates for the Ritz |
!        | values in both systems, need the last row of the    |
!        | eigenvector matrix. Remember, it`s in factored form |
!        %-----------------------------------------------------%
!
         do 65 j = 1, ncv-1
            workl(ihb+j-1) = zero
  65     continue
         workl(ihb+ncv-1) = one
         call sorm2r('Left', 'Transpose'  , ncv       ,
     &                1     , nconv        , workl(iq) ,
     &                ldq   , workl(iw+ncv), workl(ihb),
     &                ncv   , temp         , ierr      )
!
      else if (rvec .and. howmny .eq. 'S') then
!
!     Not yet implemented. See remark 2 above.
!
      end if
!
      if (type .eq. 'REGULR' .and. rvec) then
!
            do 70 j=1, ncv
               workl(ihb+j-1) = rnorm * abs( workl(ihb+j-1) )
 70         continue
!
      else if (type .ne. 'REGULR' .and. rvec) then
!
!        %-------------------------------------------------%
!        | *  Determine Ritz estimates of the theta.       |
!        |    If RVEC = .true. then compute Ritz estimates |
!        |               of the theta.                     |
!        |    If RVEC = .false. then copy Ritz estimates   |
!        |              as computed by pssaupd.            |
!        | *  Determine Ritz estimates of the lambda.      |
!        %-------------------------------------------------%
!
         call sscal (ncv, bnorm2, workl(ihb), 1)
         if (type .eq. 'SHIFTI') then
!
            do 80 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) )
     &                        / workl(iw+k-1)**2
 80         continue
!
         else if (type .eq. 'BUCKLE') then
!
            do 90 k=1, ncv
               workl(ihb+k-1) = sigma * abs( workl(ihb+k-1) )
     &                        / ( workl(iw+k-1)-one )**2
 90         continue
!
         else if (type .eq. 'CAYLEY') then
!
            do 100 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1)
     &                        / workl(iw+k-1)*(workl(iw+k-1)-one) )
 100        continue
!
         end if
!
      end if
!
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call psvout (comm, logfil, nconv, d, ndigit,
     &          '_seupd: Untransformed converged Ritz values')
         call psvout (comm, logfil, nconv, workl(ihb), ndigit,
     &     '_seupd: Ritz estimates of the untransformed Ritz values')
      else if (msglvl .gt. 1) then
         call psvout (comm, logfil, nconv, d, ndigit,
     &          '_seupd: Converged Ritz values')
         call psvout (comm, logfil, nconv, workl(ihb), ndigit,
     &     '_seupd: Associated Ritz estimates')
      end if
!
!     %-------------------------------------------------%
!     | Ritz vector purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 3,4,5. See reference 7               |
!     %-------------------------------------------------%
!
      if (rvec .and. (type .eq. 'SHIFTI' .or. type .eq. 'CAYLEY')) then
!
         do 110 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / workl(iw+k)
 110     continue
!
      else if (rvec .and. type .eq. 'BUCKLE') then
!
         do 120 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1)
     &                  / (workl(iw+k)-one)
 120     continue
!
      end if
!
      if (type .ne. 'REGULR')
     &   call sger(n, nconv, one, resid, 1, workl(iw), 1, z, ldz)
!
 9000 continue
!
      return
!
!     %----------------%
!     | End of psseupd |
!     %----------------%
!
      end
