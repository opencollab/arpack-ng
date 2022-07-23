!\BeginDoc
!
!\Name: psnorm2
!
! Message Passing Layer: MPI
!
!\Description:
!
!\Usage:
!  call psnorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
!
!-----------------------------------------------------------------------
!
      Real function psnorm2 ( comm, n, x, inc )
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
      Real
     &             x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Real
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
      Real
     &             snrm2
      External     snrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      psnorm2 = snrm2( n, x, inc)
!
      buf = psnorm2
      call MPI_ALLREDUCE( [buf], buf2, 1, MPI_REAL,
     &                    MPI_MAX, comm, ierr )
      max = buf2(1)
      if ( max .eq. zero ) then
         psnorm2 = zero
      else
         buf = (psnorm2/max)**2.0
         call MPI_ALLREDUCE( [buf], buf2, 1, MPI_REAL,
     &                       MPI_SUM, comm, ierr )
         psnorm2 = max * sqrt(abs(buf2(1)))
      endif
!
!     %----------------%
!     | End of psnorm2 |
!     %----------------%
!
      return
      end
