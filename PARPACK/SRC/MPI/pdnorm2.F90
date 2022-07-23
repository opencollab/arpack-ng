!\BeginDoc
!
!\Name: pdnorm2
!
! Message Passing Layer: MPI
!
!\Description:
!
!\Usage:
!  call pdnorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
!
!-----------------------------------------------------------------------
!
      Double precision function pdnorm2 ( comm, n, x, inc )
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
      Double precision
     &             x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Double precision
     &             max, buf, zero, buf2(1)
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
      Double precision
     &             dnrm2
      External     dnrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      pdnorm2 = dnrm2( n, x, inc)
!
      buf = pdnorm2
      call MPI_ALLREDUCE( [buf], buf2, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, comm, ierr )
      max = buf2(1)
      if ( max .eq. zero ) then
         pdnorm2 = zero
      else
         buf = (pdnorm2/max)**2.0
         call MPI_ALLREDUCE( [buf], buf2, 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, comm, ierr )
         pdnorm2 = max * sqrt(abs(buf2(1)))
      endif
!
!     %----------------%
!     | End of pdnorm2 |
!     %----------------%
!
      return
      end
