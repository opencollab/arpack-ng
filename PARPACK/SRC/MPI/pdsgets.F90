!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: pdsgets
!
! Message Passing Layer: MPI
!
!\Description:
!  Given the eigenvalues of the symmetric tridiagonal matrix H,
!  computes the NP shifts AMU that are zeros of the polynomial of
!  degree NP which filters out components of the unwanted eigenvectors
!  corresponding to the AMU's based on some given criteria.
!
!  NOTE: This is called even in the case of user specified shifts in
!  order to sort the eigenvalues, and error bounds of H for later use.
!
!\Usage:
!  call pdsgets
!     ( COMM, ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
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
!          'LM' -> KEV eigenvalues of largest magnitude are retained.
!          'SM' -> KEV eigenvalues of smallest magnitude are retained.
!          'LA' -> KEV eigenvalues of largest value are retained.
!          'SA' -> KEV eigenvalues of smallest value are retained.
!          'BE' -> KEV eigenvalues, half from each end of the spectrum.
!                  If KEV is odd, compute one more from the high end.
!
!  KEV      Integer.  (INPUT)
!          KEV+NP is the size of the matrix H.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be computed.
!
!  RITZ    Double precision array of length KEV+NP.  (INPUT/OUTPUT)
!          On INPUT, RITZ contains the eigenvalues of H.
!          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
!          are in the first NP locations and the wanted part is in
!          the last KEV locations.  When exact shifts are selected, the
!          unwanted part corresponds to the shifts to be applied.
!
!  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!  SHIFTS  Double precision array of length NP.  (INPUT/OUTPUT)
!          On INPUT:  contains the user specified shifts if ISHIFT = 0.
!          On OUTPUT: contains the shifts sorted into decreasing order
!          of magnitude with respect to the Ritz estimates contained in
!          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
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
!     dsortr  ARPACK utility sorting routine.
!     pivout  Parallel ARPACK utility routine that prints integers.
!     arscnd  ARPACK utility routine for timing.
!     pdvout  Parallel ARPACK utility routine that prints vectors.
!     dcopy   Level 1 BLAS that copies one vector to another.
!     dswap   Level 1 BLAS that swaps the contents of two vectors.
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
!     Starting Point: Serial Code FILE: sgets.F   SID: 2.3
!
!\SCCS Information:
! FILE: sgets.F   SID: 1.2   DATE OF SID: 2/22/96
!
!\Remarks
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pdsgets
     &      ( comm, ishift, which, kev, np, ritz, bounds, shifts )
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
     &           bounds(kev+np), ritz(kev+np), shifts(np)
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
      integer    kevd2, msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   dswap, dcopy, dsortr, arscnd
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    max, min
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
      msglvl = msgets
!
      if (which .eq. 'BE') then
!
!        %-----------------------------------------------------%
!        | Both ends of the spectrum are requested.            |
!        | Sort the eigenvalues into algebraically increasing  |
!        | order first then swap high end of the spectrum next |
!        | to low end in appropriate locations.                |
!        | NOTE: when np < floor(kev/2) be careful not to swap |
!        | overlapping locations.                              |
!        %-----------------------------------------------------%
!
         call dsortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2
         if ( kev .gt. 1 ) then
            call dswap ( min(kevd2,np), ritz, 1,
     &                   ritz( max(kevd2,np)+1 ), 1)
            call dswap ( min(kevd2,np), bounds, 1,
     &                   bounds( max(kevd2,np)+1 ), 1)
         end if
!
      else
!
!        %----------------------------------------------------%
!        | LM, SM, LA, SA case.                               |
!        | Sort the eigenvalues of H into the desired order   |
!        | and apply the resulting order to BOUNDS.           |
!        | The eigenvalues are sorted so that the wanted part |
!        | are always in the last KEV locations.              |
!        %----------------------------------------------------%
!
         call dsortr (which, .true., kev+np, ritz, bounds)
      end if
!
      if (ishift .eq. 1 .and. np .gt. 0) then
!
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first.       |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when the shifts  |
!        | are applied in subroutine pdsapps.                    |
!        %-------------------------------------------------------%
!
         call dsortr ('SM', .true., np, bounds, ritz)
         call dcopy (np, ritz, 1, shifts, 1)
      end if
!
      call arscnd (t1)
      tsgets = tsgets + (t1 - t0)
!
      if (msglvl .gt. 0) then
         call pivout (comm, logfil, 1, [kev], ndigit, '_sgets: KEV is')
         call pivout (comm, logfil, 1, [np], ndigit, '_sgets: NP is')
         call pdvout (comm, logfil, kev+np, ritz, ndigit,
     &        '_sgets: Eigenvalues of current H matrix')
         call pdvout (comm, logfil, kev+np, bounds, ndigit,
     &        '_sgets: Associated Ritz estimates')
      end if
!
      return
!
!     %----------------%
!     | End of pdsgets |
!     %----------------%
!
      end
