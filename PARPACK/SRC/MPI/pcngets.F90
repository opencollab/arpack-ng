!\BeginDoc
!
!\Name: pcngets
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
!  call pcngets
!      ( COMM, ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS )
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
!          'LR' -> want the KEV eigenvalues of largest REAL part.
!          'SR' -> want the KEV eigenvalues of smallest REAL part.
!          'LI' -> want the KEV eigenvalues of largest imaginary part.
!          'SI' -> want the KEV eigenvalues of smallest imaginary part.
!
!  KEV     Integer.  (INPUT)
!          The number of desired eigenvalues.
!
!  NP      Integer.  (INPUT)
!          The number of shifts to compute.
!
!  RITZ    Complex array of length KEV+NP.  (INPUT/OUTPUT)
!          On INPUT, RITZ contains the the eigenvalues of H.
!          On OUTPUT, RITZ are sorted so that the unwanted
!          eigenvalues are in the first NP locations and the wanted
!          portion is in the last KEV locations.  When exact shifts are
!          selected, the unwanted part corresponds to the shifts to
!          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
!          are further sorted so that the ones with largest Ritz values
!          are first.
!
!  BOUNDS  Complex array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  Complex
!
!\Routines called:
!     csortc  ARPACK sorting routine.
!     pivout  Parallel ARPACK utility routine that prints integers.
!     arscnd  ARPACK utility routine for timing.
!     pcvout  Parallel ARPACK utility routine that prints vectors.
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
!     Starting Point: Serial Complex Code FILE: ngets.F   SID: 2.1
!
!\SCCS Information:
! FILE: ngets.F   SID: 1.2   DATE OF SID: 4/19/96
!
!\Remarks
!     1. This routine does not keep complex conjugate pairs of
!        eigenvalues together.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pcngets ( comm, ishift, which, kev, np, ritz, bounds)
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
      Complex
     &           bounds(kev+np), ritz(kev+np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex
     &           one, zero
      parameter (one = (1.0, 0.0), zero = (0.0, 0.0))
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
      external   pcvout,  csortc, arscnd
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
      msglvl = mcgets
!
      call csortc (which, .true., kev+np, ritz, bounds)
!
      if ( ishift .eq. 1 ) then
!
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first        |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when the shifts  |
!        | are applied in subroutine pcnapps.                    |
!        | Be careful and use 'SM' since we want to sort BOUNDS! |
!        %-------------------------------------------------------%
!
         call csortc ( 'SM', .true., np, bounds, ritz )
!
      end if
!
      call arscnd (t1)
      tcgets = tcgets + (t1 - t0)
!
      if (msglvl .gt. 0) then
         call pivout (comm, logfil, 1, [kev], ndigit, '_ngets: KEV is')
         call pivout (comm, logfil, 1, [np], ndigit, '_ngets: NP is')
         call pcvout (comm, logfil, kev+np, ritz, ndigit,
     &        '_ngets: Eigenvalues of current H matrix ')
         call pcvout (comm, logfil, kev+np, bounds, ndigit,
     &      '_ngets: Ritz estimates of the current KEV+NP Ritz values')
      end if
!
      return
!
!     %----------------%
!     | End of pcngets |
!     %----------------%
!
      end
