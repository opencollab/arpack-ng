!\BeginDoc
!
!\Name: pscnorm2
!
! Message Passing Layer: MPI
!
!\Description:
!
!\Usage:
!  call pscnorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
!
!-----------------------------------------------------------------------
!
      Real function pscnorm2 ( comm, n, x, inc )
!
      include   'mpif.h'
!
!     %---------------%
!     | MPI Variables |
!     %---------------%
!
      integer    comm, ierr
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer      n, inc
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Complex
     &             x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Real
     &             max(1), buf, zero
      parameter    ( zero = 0.0 )
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &             scnrm2
      External     scnrm2
      Real         buf2(1)
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      pscnorm2 = scnrm2( n, x, inc)
!
      buf = pscnorm2
      call MPI_ALLREDUCE( [buf], max, 1, MPI_REAL,
     &                    MPI_MAX, comm, ierr )
      if ( max(1) .eq. zero ) then
         pscnorm2 = zero
      else
         buf = (pscnorm2/max(1))**2.0
         call MPI_ALLREDUCE( [buf], buf2, 1, MPI_REAL,
     &                       MPI_SUM, comm, ierr )
         pscnorm2 = max(1) * sqrt(abs(buf2(1)))
      endif
!
!     %-----------------%
!     | End of pscnorm2 |
!     %-----------------%
!
      return
      end
