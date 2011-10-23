c\BeginDoc
c
c\Name: pdnorm2
c
c Message Passing Layer: BLACS
c
c\Description:
c
c\Usage:
c  call pdnorm2 ( COMM, N, X, INC )
c
c\Arguments
c  COMM    BLACS Communicator for the processor grid.  (INPUT)
c
c\SCCS Information:
c FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
c
c-----------------------------------------------------------------------
c
      Double precision function pdnorm2 ( comm, n, x, inc )
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
      Double precision
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
     &             dnrm2
      External     dnrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      pdnorm2 = dnrm2( n, x, inc)
c
      max = pdnorm2
      call dgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         pdnorm2 = zero
      else
         pdnorm2 = (pdnorm2/max)**2.0
         call dgsum2d( comm, 'All', ' ', 1, 1, pdnorm2, 1, -1, -1 )
         pdnorm2 = max * sqrt(abs(pdnorm2))
      endif
c
c     %----------------%
c     | End of pdnorm2 |
c     %----------------%
c
      return
      end
