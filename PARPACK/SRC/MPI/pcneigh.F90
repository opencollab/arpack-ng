!\BeginDoc
!
!\Name: pcneigh
!
! Message Passing Layer: MPI
!
!\Description:
!  Compute the eigenvalues of the current upper Hessenberg matrix
!  and the corresponding Ritz estimates given the current residual norm.
!
!\Usage:
!  call pcneigh
!     ( COMM, RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  RNORM   Real scalar.  (INPUT)
!          Residual norm corresponding to the current upper Hessenberg
!          matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the matrix H.
!
!  H       Complex N by N array.  (INPUT)
!          H contains the current upper Hessenberg matrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZ    Complex array of length N.  (OUTPUT)
!          On output, RITZ(1:N) contains the eigenvalues of H.
!
!  BOUNDS  Complex array of length N.  (OUTPUT)
!          On output, BOUNDS contains the Ritz estimates associated with
!          the eigenvalues held in RITZ.  This is equal to RNORM
!          times the last components of the eigenvectors corresponding
!          to the eigenvalues in RITZ.
!
!  Q       Complex N by N array.  (WORKSPACE)
!          Workspace needed to store the eigenvectors of H.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Complex work array of length N**2 + 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  This is needed to keep the full Schur form
!          of H and also in the calculation of the eigenvectors of H.
!
!  RWORK   Real  work array of length N (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from clahqr or ctrevc.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  Complex
!
!\Routines called:
!     pivout  Parallel ARPACK utility routine that prints integers.
!     arscnd  ARPACK utility routine for timing.
!     pcmout  Parallel ARPACK utility routine that prints matrices
!     pcvout  Parallel ARPACK utility routine that prints vectors.
!     psvout  Parallel ARPACK utility routine that prints vectors.
!     clacpy  LAPACK matrix copy routine.
!     clahqr  LAPACK routine to compute the Schur form of an
!             upper Hessenberg matrix.
!     claset  LAPACK matrix initialization routine.
!     ctrevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper triangular form
!     ccopy   Level 1 BLAS that copies one vector to another.
!     csscal  Level 1 BLAS that scales a complex vector by a real number.
!     scnrm2  Level 1 BLAS that computes the norm of a vector.
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
      subroutine pcneigh (comm, rnorm, n, h, ldh, ritz, bounds,
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
      Real
     &           rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Complex
     &           bounds(n), h(ldh,n), q(ldq,n), ritz(n),
     &           workl(n*(n+3))
      Real
     &           rwork(n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex
     &           one, zero
      Real
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
      Complex
     &           vl(1)
      Real
     &           temp
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   clacpy, clahqr, csscal, ctrevc, ccopy,
     &           pcmout, pcvout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           scnrm2
      external   scnrm2
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
          call pcmout (comm, logfil, n, n, h, ldh, ndigit,
     &         '_neigh: Entering upper Hessenberg matrix H ')
      end if
!
!     %----------------------------------------------------------%
!     | 1. Compute the eigenvalues, the last components of the   |
!     |    corresponding Schur vectors and the full Schur form T |
!     |    of the current upper Hessenberg matrix H.             |
!     |    clahqr returns the full Schur form of H               |
!     |    in WORKL(1:N**2), and the Schur vectors in q.         |
!     %----------------------------------------------------------%
!
      call clacpy ('All', n, n, h, ldh, workl, n)
      call claset ('All', n, n, zero, one, q, ldq)
      call clahqr (.true., .true., n, 1, n, workl, ldh, ritz,
     &             1, n, q, ldq, ierr)
      if (ierr .ne. 0) go to 9000
!
      call ccopy (n, q(n-1,1), ldq, bounds, 1)
      if (msglvl .gt. 1) then
         call pcvout (comm, logfil, n, bounds, ndigit,
     &              '_neigh: last row of the Schur matrix for H')
      end if
!
!     %----------------------------------------------------------%
!     | 2. Compute the eigenvectors of the full Schur form T and |
!     |    apply the Schur vectors to get the corresponding      |
!     |    eigenvectors.                                         |
!     %----------------------------------------------------------%
!
      call ctrevc ('Right', 'Back', select, n, workl, n, vl, n, q,
     &             ldq, n, n, workl(n*n+1), rwork, ierr)
!
      if (ierr .ne. 0) go to 9000
!
!     %------------------------------------------------%
!     | Scale the returning eigenvectors so that their |
!     | Euclidean norms are all one. LAPACK subroutine |
!     | ctrevc returns each eigenvector normalized so  |
!     | that the element of largest magnitude has      |
!     | magnitude 1; here the magnitude of a complex   |
!     | number (x,y) is taken to be |x| + |y|.         |
!     %------------------------------------------------%
!
      do 10 j=1, n
            temp = scnrm2( n, q(1,j), 1 )
            call csscal ( n, rone / temp, q(1,j), 1 )
   10 continue
!
      if (msglvl .gt. 1) then
         call ccopy(n, q(n,1), ldq, workl, 1)
         call pcvout (comm, logfil, n, workl, ndigit,
     &              '_neigh: Last row of the eigenvector matrix for H')
      end if
!
!     %----------------------------%
!     | Compute the Ritz estimates |
!     %----------------------------%
!
      call ccopy(n, q(n,1), n, bounds, 1)
      call csscal(n, rnorm, bounds, 1)
!
      if (msglvl .gt. 2) then
         call pcvout (comm, logfil, n, ritz, ndigit,
     &              '_neigh: The eigenvalues of H')
         call pcvout (comm, logfil, n, bounds, ndigit,
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
!     | End of pcneigh |
!     %----------------%
!
      end
