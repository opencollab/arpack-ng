c\BeginDoc
c
c\Name: psnorm2
c
c Message Passing Layer: BLACS
c
c\Description:
c
c\Usage:
c  call psnorm2 ( COMM, N, X, INC )
c
c\Arguments
c  COMM    BLACS Communicator for the processor grid.  (INPUT)
c
c\SCCS Information:
c FILE: norm2.F   SID: 1.2   DATE OF SID: 2/22/96
c
c-----------------------------------------------------------------------
c
      Real function psnorm2 ( comm, n, x, inc )
c
c     %------------------------------%
c     | BLACS Variables and Routines |
c     %------------------------------%
c
      integer    comm
      external   sgsum2d, sgamx2d
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
      Real
     &             x(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      Real
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
      Real       
     &             snrm2
      External     snrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      psnorm2 = snrm2( n, x, inc)
c
      max = psnorm2
      call sgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         psnorm2 = zero
      else
         psnorm2 = (psnorm2/max)**2.0
         call sgsum2d( comm, 'All', ' ', 1, 1, psnorm2, 1, -1, -1 )
         psnorm2 = max * sqrt(abs(psnorm2))
      endif
c
c     %----------------%
c     | End of psnorm2 |
c     %----------------%
c
      return
      end
