!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: pdngets
!
! Message Passing Layer: MPI
!
!\Description:
!  Given the eigenvalues of the upper Hessenberg matrix H,
!  computes the NP shifts AMU that are zeros of the polynomial of
!  degree NP which filters out components of the unwanted eigenvectors
!  corresponding to the AMU's based on some given criteria.
!
!  NOTE: call this even in the case of user specified shifts in order
!  to sort the eigenvalues, and error bounds of H for later use.
!
!\Usage:
!  call pdngets
!     ( COMM, ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  ISHIFT  Integer.  (INPUT)
!          Method for selecting the implicit shifts at each iteration.
!          ISHIFT = 0: user specified shifts
!          ISHIFT = 1: exact shift with respect to the matrix H.
!
!  WHICH   Character*2.  (INPUT)
!          Shift selection criteria.
!          'LM' -> want the KEV eigenvalues of largest magnitude.
!          'SM' -> want the KEV eigenvalues of smallest magnitude.
!          'LR' -> want the KEV eigenvalues of largest real part.
!          'SR' -> want the KEV eigenvalues of smallest real part.
!          'LI' -> want the KEV eigenvalues of largest imaginary part.
!          'SI' -> want the KEV eigenvalues of smallest imaginary part.
!
!  KEV      Integer.  (INPUT/OUTPUT)
!           INPUT: KEV+NP is the size of the matrix H.
!           OUTPUT: Possibly increases KEV by one to keep complex conjugate
!           pairs together.
!
!  NP       Integer.  (INPUT/OUTPUT)
!           Number of implicit shifts to be computed.
!           OUTPUT: Possibly decreases NP by one to keep complex conjugate
!           pairs together.
!
!  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
!  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary
!          parts of the eigenvalues of H.
!          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
!          eigenvalues are in the first NP locations and the wanted
!          portion is in the last KEV locations.  When exact shifts are
!          selected, the unwanted part corresponds to the shifts to
!          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
!          are further sorted so that the ones with largest Ritz values
!          are first.
!
!  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
!
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     dsortc  ARPACK sorting routine.
!     dcopy   Level 1 BLAS that copies one vector to another .
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
!     Starting Point: Serial Code FILE: ngets.F   SID: 2.2
!
!\SCCS Information:
! FILE: ngets.F   SID: 1.2   DATE OF SID: 2/22/96
!
!\Remarks
!     1. xxxx
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pdngets
     &                 ( comm, ishift, which, kev, np, ritzr, ritzi,
     &                   bounds, shiftr, shifti )
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
      character*2 which
      integer    ishift, kev, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Double precision
     &           bounds(kev+np), ritzr(kev+np), ritzi(kev+np),
     &           shiftr(1), shifti(1)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   dcopy, dsortc, arscnd
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call arscnd (t0)
      msglvl = mngets
!
!     %----------------------------------------------------%
!     | LM, SM, LR, SR, LI, SI case.                       |
!     | Sort the eigenvalues of H into the desired order   |
!     | and apply the resulting order to BOUNDS.           |
!     | The eigenvalues are sorted so that the wanted part |
!     | are always in the last KEV locations.              |
!     | We first do a pre-processing sort in order to keep |
!     | complex conjugate pairs together                   |
!     %----------------------------------------------------%
!
      if (which .eq. 'LM') then
         call dsortc ('LR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SM') then
         call dsortc ('SR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LR') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SR') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LI') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SI') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      end if
!
      call dsortc (which, .true., kev+np, ritzr, ritzi, bounds)
!
!     %-------------------------------------------------------%
!     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
!     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
!     | Accordingly decrease NP by one. In other words keep   |
!     | complex conjugate pairs together.                     |
!     %-------------------------------------------------------%
!
      if (       ( ritzr(np+1) - ritzr(np) ) .eq. zero
     &     .and. ( ritzi(np+1) + ritzi(np) ) .eq. zero ) then
         np = np - 1
         kev = kev + 1
      end if
!
      if ( ishift .eq. 1 ) then
!
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first        |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when they shifts |
!        | are applied in subroutine pdnapps.                    |
!        | Be careful and use 'SR' since we want to sort BOUNDS! |
!        %-------------------------------------------------------%
!
         call dsortc ( 'SR', .true., np, bounds, ritzr, ritzi )
      end if
!
      call arscnd (t1)
      tngets = tngets + (t1 - t0)
!
      if (msglvl .gt. 0) then
         call pivout (comm, logfil, 1, [kev], ndigit, '_ngets: KEV is')
         call pivout (comm, logfil, 1, [np], ndigit, '_ngets: NP is')
         call pdvout (comm, logfil, kev+np, ritzr, ndigit,
     &        '_ngets: Eigenvalues of current H matrix -- real part')
         call pdvout (comm, logfil, kev+np, ritzi, ndigit,
     &        '_ngets: Eigenvalues of current H matrix -- imag part')
         call pdvout (comm, logfil, kev+np, bounds, ndigit,
     &      '_ngets: Ritz estimates of the current KEV+NP Ritz values')
      end if
!
      return
!
!     %----------------%
!     | End of pdngets |
!     %----------------%
!
      end
