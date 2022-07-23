!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: pdseigt
!
! Message Passing Layer: MPI
!
!\Description:
!  Compute the eigenvalues of the current symmetric tridiagonal matrix
!  and the corresponding error bounds given the current residual norm.
!
!\Usage:
!  call pdseigt
!     ( COMM, RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  RNORM   Double precision scalar.  (INPUT)
!          RNORM contains the residual norm corresponding to the current
!          symmetric tridiagonal matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the symmetric tridiagonal matrix H.
!
!  H       Double precision N by 2 array.  (INPUT)
!          H contains the symmetric tridiagonal matrix with the
!          subdiagonal in the first column starting at H(2,1) and the
!          main diagonal in second column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  EIG     Double precision array of length N.  (OUTPUT)
!          On output, EIG contains the N eigenvalues of H possibly
!          unsorted.  The BOUNDS arrays are returned in the
!          same sorted order as EIG.
!
!  BOUNDS  Double precision array of length N.  (OUTPUT)
!          On output, BOUNDS contains the error estimates corresponding
!          to the eigenvalues EIG.  This is equal to RNORM times the
!          last components of the eigenvectors corresponding to the
!          eigenvalues in EIG.
!
!  WORKL   Double precision work array of length 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from dstqrb.
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
!     dstqrb  ARPACK routine that computes the eigenvalues and the
!             last components of the eigenvectors of a symmetric
!             and tridiagonal matrix.
!     arscnd  ARPACK utility routine for timing.
!     pdvout  Parallel ARPACK utility routine that prints vectors.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     dscal   Level 1 BLAS that scales a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Code FILE: seigt.F   SID: 2.2
!
!\SCCS Information:
! FILE: seigt.F   SID: 1.3   DATE OF SID: 4/19/96
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pdseigt
     &   ( comm, rnorm, n, h, ldh, eig, bounds, workl, ierr )
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
      integer    ierr, ldh, n
      Double precision
     &           rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Double precision
     &           eig(n), bounds(n), h(ldh,2), workl(3*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &           zero
      parameter (zero = 0.0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, k, msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   dcopy, dstqrb, pdvout, arscnd
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
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call arscnd (t0)
      msglvl = mseigt
!
      if (msglvl .gt. 0) then
         call pdvout (comm, logfil, n, h(1,2), ndigit,
     &              '_seigt: main diagonal of matrix H')
         if (n .gt. 1) then
         call pdvout (comm, logfil, n-1, h(2,1), ndigit,
     &              '_seigt: sub diagonal of matrix H')
         end if
      end if
!
!
      call dcopy  (n, h(1,2), 1, eig, 1)
      call dcopy  (n-1, h(2,1), 1, workl, 1)
      call dstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr .ne. 0) go to 9000
!
      if (msglvl .gt. 1) then
         call pdvout (comm, logfil, n, bounds, ndigit,
     &              '_seigt: last row of the eigenvector matrix for H')
      end if
!
!     %-----------------------------------------------%
!     | Finally determine the error bounds associated |
!     | with the n Ritz values of H.                  |
!     %-----------------------------------------------%
!
      do 30 k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
   30 continue
!
      call arscnd (t1)
      tseigt = tseigt + (t1 - t0)
!
 9000 continue
      return
!
!     %-----------------%
!     |  End of pdseigt |
!     %-----------------%
!
      end
