!\BeginDoc
!
!\Name: pdznorm2
!
! Message Passing Layer: MPI
!
!\Description:
!
!\Usage:
!  call pdznorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
!
!-----------------------------------------------------------------------
!
      Double precision function pdznorm2 ( comm, n, x, inc )
!
#ifdef HAVE_MPI_ICB
      use :: mpi_f08
#else
#include "mpif.h"
#endif

!
!     %---------------%
!     | MPI Variables |
!     %---------------%
!
#ifdef HAVE_MPI_ICB
      type(MPI_Comm) comm
      integer*4 ierr
#else
      integer    comm, ierr
#endif
      integer*4  cnt

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
      Complex*16&
                   x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Double precision&
                   max(1), buf, zero
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
      Double precision&
                   dznrm2, buf2(1)
      External     dznrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      pdznorm2 = dznrm2( n, x, inc)
!
      buf = pdznorm2
      cnt = 1
      call MPI_ALLREDUCE( [buf], max, cnt, MPI_DOUBLE_PRECISION,&
                          MPI_MAX, comm, ierr )
      if ( max(1) .eq. zero ) then
         pdznorm2 = zero
      else
         buf = (pdznorm2/max(1))**2.0
         cnt = 1
         call MPI_ALLREDUCE( [buf], buf2, cnt, MPI_DOUBLE_PRECISION,&
                             MPI_SUM, comm, ierr )
         pdznorm2 = max(1) * sqrt(abs(buf2(1)))
      endif
!
!     %-----------------%
!     | End of pdznorm2 |
!     %-----------------%
!
      return
      end
