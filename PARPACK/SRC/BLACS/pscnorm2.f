c\BeginDoc
c
c\Name: pscnorm2
c
c Message Passing Layer: BLACS
c
c\Description:
c
c\Usage:
c  call pscnorm2 ( COMM, N, X, INC )
c
c\Arguments
c  COMM    BLACS Communicator for the processor grid.  (INPUT)
c
c\SCCS Information:
c FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
c
c-----------------------------------------------------------------------
c
      Real function pscnorm2 ( comm, n, x, inc )
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
      Complex
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
     &             scnrm2
      External     scnrm2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      pscnorm2 = scnrm2( n, x, inc)
c
      max = pscnorm2
      call sgamx2d( comm, 'All', ' ', 1, 1, max, 1, ra, ca,
     &              -1, -1, -1 )
      if ( max .eq. zero ) then
         pscnorm2 = zero
      else
         pscnorm2 = (pscnorm2/max)**2.0
         call sgsum2d( comm, 'All', ' ', 1, 1, pscnorm2, 1, -1, -1 )
         pscnorm2 = max * sqrt(abs(pscnorm2))
      endif
c
c     %-----------------%
c     | End of pscnorm2 |
c     %-----------------%
c
      return
      end
