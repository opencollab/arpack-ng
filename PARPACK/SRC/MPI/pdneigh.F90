!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: pdneigh
!
! Message Passing Layer: MPI
!
!\Description:
!  Compute the eigenvalues of the current upper Hessenberg matrix
!  and the corresponding Ritz estimates given the current residual norm.
!
!\Usage:
!  call pdneigh
!     ( COMM, RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  RNORM   Double precision scalar.  (INPUT)
!          Residual norm corresponding to the current upper Hessenberg
!          matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the matrix H.
!
!  H       Double precision N by N array.  (INPUT)
!          H contains the current upper Hessenberg matrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZR,  Double precision arrays of length N.  (OUTPUT)
!  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real
!          (respectively imaginary) parts of the eigenvalues of H.
!
!  BOUNDS  Double precision array of length N.  (OUTPUT)
!          On output, BOUNDS contains the Ritz estimates associated with
!          the eigenvalues RITZR and RITZI.  This is equal to RNORM
!          times the last components of the eigenvectors corresponding
!          to the eigenvalues in RITZR and RITZI.
!
!  Q       Double precision N by N array.  (WORKSPACE)
!          Workspace needed to store the eigenvectors of H.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  This is needed to keep the full Schur form
!          of H and also in the calculation of the eigenvectors of H.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from dlaqrb or dtrevc.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     dlahqr  ARPACK routine to compute the real Schur form of an
!             upper Hessenberg matrix and last row of the Schur vectors.
!     arscnd  ARPACK utility routine for timing.
!     dmout   ARPACK utility routine that prints matrices
!     dvout   ARPACK utility routine that prints vectors.
!     dlacpy  LAPACK matrix copy routine.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper quasi-triangular form
!     dgemv   Level 2 BLAS routine for matrix vector multiplication.
!     dcopy   Level 1 BLAS that copies one vector to another .
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     dscal   Level 1 BLAS that scales a vector.
!
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
!     Starting Point: Serial Code FILE: neigh.F   SID: 2.2
!
!\SCCS Information:
! FILE: neigh.F   SID: 1.2   DATE OF SID: 2/22/96
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pdneigh ( comm, rnorm, n, h, ldh, ritzr, ritzi, bounds,
     &                    q, ldq, workl, ierr)
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
      integer    ierr, n, ldh, ldq
      Double precision
     &           rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Double precision
     &           bounds(n), h(ldh,n), q(ldq,n), ritzi(n), ritzr(n),
     &           workl(n*(n+3))
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    select(1)
      integer    i, iconj, msglvl
      Double precision
     &           temp, vl(1)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   dcopy, dlacpy, dlahqr, dtrevc, dvout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &           dlapy2, dnrm2
      external   dlapy2, dnrm2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic  abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call arscnd (t0)
      msglvl = mneigh
!
      if (msglvl .gt. 2) then
          call pdmout (comm, logfil, n, n, h, ldh, ndigit,
     &         '_neigh: Entering upper Hessenberg matrix H ')
      end if
!
!     %-----------------------------------------------------------%
!     | 1. Compute the eigenvalues, the last components of the    |
!     |    corresponding Schur vectors and the full Schur form T  |
!     |    of the current upper Hessenberg matrix H.              |
!     | dlahqr returns the full Schur form of H in WORKL(1:N**2)  |
!     | and the last components of the Schur vectors in BOUNDS.   |
!     %-----------------------------------------------------------%
!
      call dlacpy ('All', n, n, h, ldh, workl, n)
      do 5 j = 1, n-1
          bounds(j) = zero
   5  continue
      bounds(n) = 1
      call dlahqr(.true., .true., n, 1, n, workl, n, ritzr, ritzi, 1, 1,
     &            bounds, 1, ierr)
      if (ierr .ne. 0) go to 9000
!
      if (msglvl .gt. 1) then
         call pdvout (comm, logfil, n, bounds, ndigit,
     &              '_neigh: last row of the Schur matrix for H')
      end if
!
!     %-----------------------------------------------------------%
!     | 2. Compute the eigenvectors of the full Schur form T and  |
!     |    apply the last components of the Schur vectors to get  |
!     |    the last components of the corresponding eigenvectors. |
!     | Remember that if the i-th and (i+1)-st eigenvalues are    |
!     | complex conjugate pairs, then the real & imaginary part   |
!     | of the eigenvector components are split across adjacent   |
!     | columns of Q.                                             |
!     %-----------------------------------------------------------%
!
      call dtrevc ('R', 'A', select, n, workl, n, vl, n, q, ldq,
     &             n, n, workl(n*n+1), ierr)
!
      if (ierr .ne. 0) go to 9000
!
!     %------------------------------------------------%
!     | Scale the returning eigenvectors so that their |
!     | euclidean norms are all one. LAPACK subroutine |
!     | dtrevc returns each eigenvector normalized so  |
!     | that the element of largest magnitude has      |
!     | magnitude 1; here the magnitude of a complex   |
!     | number (x,y) is taken to be |x| + |y|.         |
!     %------------------------------------------------%
!
      iconj = 0
      do 10 i=1, n
         if ( abs( ritzi(i) ) .le. zero ) then
!
!           %----------------------%
!           | Real eigenvalue case |
!           %----------------------%
!
            temp = dnrm2( n, q(1,i), 1 )
            call dscal ( n, one / temp, q(1,i), 1 )
         else
!
!           %-------------------------------------------%
!           | Complex conjugate pair case. Note that    |
!           | since the real and imaginary part of      |
!           | the eigenvector are stored in consecutive |
!           | columns, we further normalize by the      |
!           | square root of two.                       |
!           %-------------------------------------------%
!
            if (iconj .eq. 0) then
               temp = dlapy2( dnrm2( n, q(1,i), 1 ),
     &                        dnrm2( n, q(1,i+1), 1 ) )
               call dscal ( n, one / temp, q(1,i), 1 )
               call dscal ( n, one / temp, q(1,i+1), 1 )
               iconj = 1
            else
               iconj = 0
            end if
         end if
   10 continue
!
      call dgemv ('T', n, n, one, q, ldq, bounds, 1, zero, workl, 1)
!
      if (msglvl .gt. 1) then
         call pdvout (comm, logfil, n, workl, ndigit,
     &              '_neigh: Last row of the eigenvector matrix for H')
      end if
!
!     %----------------------------%
!     | Compute the Ritz estimates |
!     %----------------------------%
!
      iconj = 0
      do 20 i = 1, n
         if ( abs( ritzi(i) ) .le. zero ) then
!
!           %----------------------%
!           | Real eigenvalue case |
!           %----------------------%
!
            bounds(i) = rnorm * abs( workl(i) )
         else
!
!           %-------------------------------------------%
!           | Complex conjugate pair case. Note that    |
!           | since the real and imaginary part of      |
!           | the eigenvector are stored in consecutive |
!           | columns, we need to take the magnitude    |
!           | of the last components of the two vectors |
!           %-------------------------------------------%
!
            if (iconj .eq. 0) then
               bounds(i) = rnorm * dlapy2( workl(i), workl(i+1) )
               bounds(i+1) = bounds(i)
               iconj = 1
            else
               iconj = 0
            end if
         end if
   20 continue
!
      if (msglvl .gt. 2) then
         call pdvout (comm, logfil, n, ritzr, ndigit,
     &              '_neigh: Real part of the eigenvalues of H')
         call pdvout (comm, logfil, n, ritzi, ndigit,
     &              '_neigh: Imaginary part of the eigenvalues of H')
         call pdvout (comm, logfil, n, bounds, ndigit,
     &              '_neigh: Ritz estimates for the eigenvalues of H')
      end if
!
      call arscnd (t1)
      tneigh = tneigh + (t1 - t0)
!
 9000 continue
      return
!
!     %----------------%
!     | End of pdneigh |
!     %----------------%
!
      end
