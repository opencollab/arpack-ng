!\BeginDoc
!
!\Name: pzneigh
!
! Message Passing Layer: MPI
!
!\Description:
!  Compute the eigenvalues of the current upper Hessenberg matrix
!  and the corresponding Ritz estimates given the current residual norm.
!
!\Usage:
!  call pzneigh
!     ( COMM, RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
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
!  H       Complex*16 N by N array.  (INPUT)
!          H contains the current upper Hessenberg matrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZ    Complex*16 array of length N.  (OUTPUT)
!          On output, RITZ(1:N) contains the eigenvalues of H.
!
!  BOUNDS  Complex*16 array of length N.  (OUTPUT)
!          On output, BOUNDS contains the Ritz estimates associated with
!          the eigenvalues held in RITZ.  This is equal to RNORM
!          times the last components of the eigenvectors corresponding
!          to the eigenvalues in RITZ.
!
!  Q       Complex*16 N by N array.  (WORKSPACE)
!          Workspace needed to store the eigenvectors of H.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Complex*16 work array of length N**2 + 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  This is needed to keep the full Schur form
!          of H and also in the calculation of the eigenvectors of H.
!
!  RWORK   Double precision  work array of length N (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from zlahqr or ztrevc.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  Complex*16
!
!\Routines called:
!     pivout  Parallel ARPACK utility routine that prints integers.
!     arscnd  ARPACK utility routine for timing.
!     pzmout  Parallel ARPACK utility routine that prints matrices
!     pzvout  Parallel ARPACK utility routine that prints vectors.
!     pdvout  Parallel ARPACK utility routine that prints vectors.
!     zlacpy  LAPACK matrix copy routine.
!     zlahqr  LAPACK routine to compute the Schur form of an
!             upper Hessenberg matrix.
!     zlaset  LAPACK matrix initialization routine.
!     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper triangular form
!     zcopy   Level 1 BLAS that copies one vector to another.
!     zdscal  Level 1 BLAS that scales a complex vector by a real number.
!     dznrm2  Level 1 BLAS that computes the norm of a vector.
!
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
!     Starting Point: Serial Complex Code FILE: neigh.F   SID: 2.1
!
!\SCCS Information:
! FILE: neigh.F   SID: 1.2   DATE OF SID: 4/19/96
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pzneigh (comm, rnorm, n, h, ldh, ritz, bounds,
     &                   q, ldq, workl, rwork, ierr)
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
      Complex*16
     &           bounds(n), h(ldh,n), q(ldq,n), ritz(n),
     &           workl(n*(n+3))
      Double precision
     &           rwork(n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex*16
     &           one, zero
      Double precision
     &           rone
      parameter  (one = (1.0, 0.0), zero = (0.0, 0.0),
     &           rone = 1.0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    select(1)
      integer    j,  msglvl
      Complex*16
     &           vl(1)
      Double precision
     &           temp
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   zlacpy, zlahqr, zdscal, ztrevc, zcopy,
     &           pzmout, pzvout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &           dznrm2
      external   dznrm2
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
      msglvl = mceigh
!
      if (msglvl .gt. 2) then
          call pzmout (comm, logfil, n, n, h, ldh, ndigit,
     &         '_neigh: Entering upper Hessenberg matrix H ')
      end if
!
!     %----------------------------------------------------------%
!     | 1. Compute the eigenvalues, the last components of the   |
!     |    corresponding Schur vectors and the full Schur form T |
!     |    of the current upper Hessenberg matrix H.             |
!     |    zlahqr returns the full Schur form of H               |
!     |    in WORKL(1:N**2), and the Schur vectors in q.         |
!     %----------------------------------------------------------%
!
      call zlacpy ('All', n, n, h, ldh, workl, n)
      call zlaset ('All', n, n, zero, one, q, ldq)
      call zlahqr (.true., .true., n, 1, n, workl, ldh, ritz,
     &             1, n, q, ldq, ierr)
      if (ierr .ne. 0) go to 9000
!
      call zcopy (n, q(n-1,1), ldq, bounds, 1)
      if (msglvl .gt. 1) then
         call pzvout (comm, logfil, n, bounds, ndigit,
     &              '_neigh: last row of the Schur matrix for H')
      end if
!
!     %----------------------------------------------------------%
!     | 2. Compute the eigenvectors of the full Schur form T and |
!     |    apply the Schur vectors to get the corresponding      |
!     |    eigenvectors.                                         |
!     %----------------------------------------------------------%
!
      call ztrevc ('Right', 'Back', select, n, workl, n, vl, n, q,
     &             ldq, n, n, workl(n*n+1), rwork, ierr)
!
      if (ierr .ne. 0) go to 9000
!
!     %------------------------------------------------%
!     | Scale the returning eigenvectors so that their |
!     | Euclidean norms are all one. LAPACK subroutine |
!     | ztrevc returns each eigenvector normalized so  |
!     | that the element of largest magnitude has      |
!     | magnitude 1; here the magnitude of a complex   |
!     | number (x,y) is taken to be |x| + |y|.         |
!     %------------------------------------------------%
!
      do 10 j=1, n
            temp = dznrm2( n, q(1,j), 1 )
            call zdscal ( n, rone / temp, q(1,j), 1 )
   10 continue
!
      if (msglvl .gt. 1) then
         call zcopy(n, q(n,1), ldq, workl, 1)
         call pzvout (comm, logfil, n, workl, ndigit,
     &              '_neigh: Last row of the eigenvector matrix for H')
      end if
!
!     %----------------------------%
!     | Compute the Ritz estimates |
!     %----------------------------%
!
      call zcopy(n, q(n,1), n, bounds, 1)
      call zdscal(n, rnorm, bounds, 1)
!
      if (msglvl .gt. 2) then
         call pzvout (comm, logfil, n, ritz, ndigit,
     &              '_neigh: The eigenvalues of H')
         call pzvout (comm, logfil, n, bounds, ndigit,
     &              '_neigh: Ritz estimates for the eigenvalues of H')
      end if
!
      call arscnd(t1)
      tceigh = tceigh + (t1 - t0)
!
 9000 continue
      return
!
!     %----------------%
!     | End of pzneigh |
!     %----------------%
!
      end
