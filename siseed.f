c     Set iseed.
      subroutine siseed(seed1_c, seed2_c, seed3_c, seed4_c)
c
      integer seed1_c, seed2_c, seed3_c, seed4_c
c
      integer iseed(4)
      logical inits
c     Save: preserves items in a subprogram after the RETURN
      save iseed, inits
c
c     Set iseed.
      iseed(1) = seed1_c
      iseed(2) = seed2_c
      iseed(3) = seed3_c
      iseed(4) = seed4_c
c     We want to use these seeds: no need to initialise.
      inits = .false.
c     Save: preserves items in a subprogram after the RETURN
      return
      end
