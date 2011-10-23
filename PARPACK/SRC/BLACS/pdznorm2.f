c\BeginDoc
c
c\Name: pdznorm2
c
c Message Passing Layer: BLACS
c
c\Description:
c
c\Usage:
c  call pdznorm2 ( COMM, N, X, INC )
c
c\Arguments
c  COMM    BLACS Communicator for the processor grid.  (INPUT)
c
c\SCCS Information:
c FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
c
c-----------------------------------------------------------------------
c
      Double precision function pdznorm2 ( comm, n, x, inc )
c
c     %------------------------------%
c     | BLACS Variables and Routines |
c     %------------------------------%
c
      integer    comm
      external   dgsum2d, dgamx2d
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer      n, inc
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16
     &             x(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      Double precision
     &             max, buf, zero
      parameter    ( zero = 0.0 )
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision       
     &             dznrm2
      External     dznrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      pdznorm2 = dznrm2( n, x, inc)
c
      max = pdznorm2
      call dgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         pdznorm2 = zero
      else
         pdznorm2 = (pdznorm2/max)**2.0
         call dgsum2d( comm, 'All', ' ', 1, 1, pdznorm2, 1, -1, -1 )
         pdznorm2 = max * sqrt(abs(pdznorm2))
      endif
c
c     %-----------------%
c     | End of pdznorm2 |
c     %-----------------%
c
      return
      end
