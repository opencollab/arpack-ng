!\BeginDoc
!
!\Name: pdnaup2
!
! Message Passing Layer: MPI
!
!\Description:
!  Intermediate level interface called by pdnaupd .
!
!\Usage:
!  call pdnaup2
!     ( COMM, IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
!       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
!
!\Arguments
!
!  COMM, IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in pdnaupd .
!  MODE, ISHIFT, MXITER: see the definition of IPARAM in pdnaupd .
!
!  NP      Integer.  (INPUT/OUTPUT)
!          Contains the number of implicit shifts to apply during
!          each Arnoldi iteration.
!          If ISHIFT=1, NP is adjusted dynamically at each iteration
!          to accelerate convergence and prevent stagnation.
!          This is also roughly equal to the number of matrix-vector
!          products (involving the operator OP) per Arnoldi iteration.
!          The logic for adjusting is contained within the current
!          subroutine.
!          If ISHIFT=0, NP is the number of shifts the user needs
!          to provide via reverse communication. 0 < NP < NCV-NEV.
!          NP may be less than NCV-NEV for two reasons. The first, is
!          to keep complex conjugate pairs of "wanted" Ritz values
!          together. The second, is that a leading block of the current
!          upper Hessenberg matrix has split off and contains "unwanted"
!          Ritz values.
!          Upon termination of the IRA iteration, NP contains the number
!          of "converged" wanted Ritz values.
!
!  IUPD    Integer.  (INPUT)
!          IUPD .EQ. 0: use explicit restart instead implicit update.
!          IUPD .NE. 0: use implicit update.
!
!  V       Double precision  N by (NEV+NP) array.  (INPUT/OUTPUT)
!          The Arnoldi basis vectors are returned in the first NEV
!          columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Double precision  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
!          H is used to store the generated upper Hessenberg matrix
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZR,  Double precision  arrays of length NEV+NP.  (OUTPUT)
!  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
!          imaginary) part of the computed Ritz values of OP.
!
!  BOUNDS  Double precision  array of length NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) contain the error bounds corresponding to
!          the computed Ritz values.
!
!  Q       Double precision  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
!          Private (replicated) work array used to accumulate the
!          rotation in the shift application step.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Double precision  work array of length at least
!          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  It is used in shifts calculation, shifts
!          application and convergence checking.
!
!          On exit, the last 3*(NEV+NP) locations of WORKL contain
!          the Ritz values (real,imaginary) and associated Ritz
!          estimates of the current Hessenberg matrix.  They are
!          listed in the same order as returned from dneigh .
!
!          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
!          of WORKL are used in reverse communication to hold the user
!          supplied shifts.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD for
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!
!  WORKD   Double precision  work array of length 3*N.  (WORKSPACE)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration !!!!!!!!!!
!          See Data Distribution Note in DNAUPD.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =     0: Normal return.
!          =     1: Maximum number of iterations taken.
!                   All possible eigenvalues of OP has been found.
!                   NP returns the number of converged Ritz values.
!          =     2: No shifts could be applied.
!          =    -8: Error return from LAPACK eigenvalue calculation;
!                   This should never happen.
!          =    -9: Starting vector is zero.
!          = -9999: Could not build an Arnoldi factorization.
!                   Size that was built in returned in NP.
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
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     pdgetv0   Parallel ARPACK initial vector generation routine.
!     pdnaitr   Parallel ARPACK Arnoldi factorization routine.
!     pdnapps   Parallel ARPACK application of implicit shifts routine.
!     dnconv    ARPACK convergence of Ritz values routine.
!     pdneigh   Parallel ARPACK compute Ritz values and error bounds routine.
!     pdngets   Parallel ARPACK reorder Ritz values and error bounds routine.
!     dsortc    ARPACK sorting routine.
!     pivout   Parallel ARPACK utility routine that prints integers.
!     arscnd   ARPACK utility routine for timing.
!     pdmout    Parallel ARPACK utility routine that prints matrices
!     pdvout    ARPACK utility routine that prints vectors.
!     pdlamch10   ScaLAPACK routine that determines machine constants.
!     dlapy2    LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dcopy     Level 1 BLAS that copies one vector to another .
!     ddot      Level 1 BLAS that computes the scalar product of two vectors.
!     pdnorm2   Parallel version of Level 1 BLAS that computes the norm of a vector.
!     dswap     Level 1 BLAS that swaps two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              Cray Research, Inc. &
!     Dept. of Computational &     CRPC / Rice University
!     Applied Mathematics          Houston, Texas
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     Starting Point: Serial Code FILE: naup2.F   SID: 2.2
!
!\SCCS Information:
! FILE: naup2.F   SID: 1.5   DATE OF SID: 06/01/00
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pdnaup2
     &   ( comm, ido, bmat, n, which, nev, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds,
     &     q, ldq, workl, ipntr, workd, info )
!
      include   'mpif.h'
!
!     %---------------%
!     | MPI Variables |
!     %---------------%
!
      integer    comm
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
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision
     &           tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(13)
      Double precision
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n),
     &           ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np),
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision
     &           one, zero
      parameter (one = 1.0 , zero = 0.0 )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  wprime*2
      logical    cnorm , getv0, initv , update, ushift
      integer    ierr  , iter , kplusp, msglvl, nconv,
     &           nevbef, nev0 , np0   , nptemp, numcnv,
     &           j
      Double precision
     &           rnorm , temp , eps23, buf2(1)
      save       cnorm , getv0, initv , update, ushift,
     &           rnorm , iter , kplusp, msglvl, nconv,
     &           nevbef, nev0 , np0   , eps23 , numcnv
!

      Double precision
     &           rnorm_buf
!
!     %-----------------------%
!     | Local array arguments |
!     %-----------------------%
!
      integer    kp(4)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   dcopy , pdgetv0 , pdnaitr , dnconv ,
     &           pdneigh , pdngets , pdnapps ,
     &           pdvout , pivout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &           ddot , pdnorm2 , dlapy2 , pdlamch10
      external   ddot , pdnorm2 , dlapy2 , pdlamch10
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min, max, abs, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido .eq. 0) then
!
         call arscnd (t0)
!
         msglvl = mnaup2
!
!        %-------------------------------------%
!        | Get the machine dependent constant. |
!        %-------------------------------------%
!
         eps23 = pdlamch10 (comm, 'Epsilon-Machine')
         eps23 = eps23**(2.0  / 3.0 )
!
         nev0   = nev
         np0    = np
!
!        %-------------------------------------%
!        | kplusp is the bound on the largest  |
!        |        Lanczos factorization built. |
!        | nconv is the current number of      |
!        |        "converged" eigenvlues.      |
!        | iter is the counter on the current  |
!        |      iteration step.                |
!        %-------------------------------------%
!
         kplusp = nev + np
         nconv  = 0
         iter   = 0
!
!        %---------------------------------------%
!        | Set flags for computing the first NEV |
!        | steps of the Arnoldi factorization.   |
!        %---------------------------------------%
!
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
!
         if (info .ne. 0) then
!
!           %--------------------------------------------%
!           | User provides the initial residual vector. |
!           %--------------------------------------------%
!
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
!
!     %---------------------------------------------%
!     | Get a possibly random starting vector and   |
!     | force it into the range of the operator OP. |
!     %---------------------------------------------%
!
   10 continue
!
      if (getv0) then
         call pdgetv0  (comm, ido, bmat, 1, initv, n, 1, v, ldv,
     &                resid, rnorm, ipntr, workd, workl, info)
!
         if (ido .ne. 99) go to 9000
!
         if (rnorm .eq. zero) then
!
!           %-----------------------------------------%
!           | The initial vector is zero. Error exit. |
!           %-----------------------------------------%
!
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
!
!     %-----------------------------------%
!     | Back from reverse communication : |
!     | continue with update step         |
!     %-----------------------------------%
!
      if (update) go to 20
!
!     %-------------------------------------------%
!     | Back from computing user specified shifts |
!     %-------------------------------------------%
!
      if (ushift) go to 50
!
!     %-------------------------------------%
!     | Back from computing residual norm   |
!     | at the end of the current iteration |
!     %-------------------------------------%
!
      if (cnorm)  go to 100
!
!     %----------------------------------------------------------%
!     | Compute the first NEV steps of the Arnoldi factorization |
!     %----------------------------------------------------------%
!
      call pdnaitr  (comm, ido, bmat, n, 0, nev, mode,
     &             resid, rnorm, v, ldv, h, ldh, ipntr,
     &             workd, workl, info)
!
!     %---------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication  |
!     | to compute operations involving OP and possibly B |
!     %---------------------------------------------------%
!
      if (ido .ne. 99) go to 9000
!
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
!     |           Each iteration implicitly restarts the Arnoldi     |
!     |           factorization in place.                            |
!     |                                                              |
!     %--------------------------------------------------------------%
!
 1000 continue
!
         iter = iter + 1
!
         if (msglvl .gt. 0) then
            call pivout (comm, logfil, 1, [iter], ndigit,
     &           '_naup2: **** Start of major iteration number ****')
         end if
!
!        %-----------------------------------------------------------%
!        | Compute NP additional steps of the Arnoldi factorization. |
!        | Adjust NP since NEV might have been updated by last call  |
!        | to the shift application routine pdnapps .                 |
!        %-----------------------------------------------------------%
!
         np  = kplusp - nev
!
         if (msglvl .gt. 1) then
            call pivout (comm, logfil, 1, [nev], ndigit,
     &     '_naup2: The length of the current Arnoldi factorization')
            call pivout (comm, logfil, 1, [np], ndigit,
     &           '_naup2: Extend the Arnoldi factorization by')
         end if
!
!        %-----------------------------------------------------------%
!        | Compute NP additional steps of the Arnoldi factorization. |
!        %-----------------------------------------------------------%
!
         ido = 0
   20    continue
         update = .true.
!
         call pdnaitr  (comm, ido, bmat, n, nev, np, mode,
     &                resid, rnorm, v, ldv,
     &                h, ldh, ipntr, workd, workl, info)
!
!        %---------------------------------------------------%
!        | ido .ne. 99 implies use of reverse communication  |
!        | to compute operations involving OP and possibly B |
!        %---------------------------------------------------%
!
         if (ido .ne. 99) go to 9000
!
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
!
         if (msglvl .gt. 1) then
            call pdvout  (comm, logfil, 1, [rnorm], ndigit,
     &           '_naup2: Corresponding B-norm of the residual')
         end if
!
!        %--------------------------------------------------------%
!        | Compute the eigenvalues and corresponding error bounds |
!        | of the current upper Hessenberg matrix.                |
!        %--------------------------------------------------------%
!
         call pdneigh  ( comm, rnorm, kplusp, h, ldh, ritzr, ritzi,
     &                  bounds, q, ldq, workl, ierr)
!
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
!
!        %----------------------------------------------------%
!        | Make a copy of eigenvalues and corresponding error |
!        | bounds obtained from pdneigh .                      |
!        %----------------------------------------------------%
!
         call dcopy (kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call dcopy (kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call dcopy (kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
!
!        %---------------------------------------------------%
!        | Select the wanted Ritz values and their bounds    |
!        | to be used in the convergence test.               |
!        | The wanted part of the spectrum and corresponding |
!        | error bounds are in the last NEV loc. of RITZR,   |
!        | RITZI and BOUNDS respectively. The variables NEV  |
!        | and NP may be updated if the NEV-th wanted Ritz   |
!        | value has a non zero imaginary part. In this case |
!        | NEV is increased by one and NP decreased by one.  |
!        | NOTE: The last two arguments of pdngets  are no    |
!        | longer used as of version 2.1.                    |
!        %---------------------------------------------------%
!
         nev = nev0
         np = np0
         numcnv = nev
         call pdngets  ( comm, ishift, which, nev, np, ritzr, ritzi,
     &                  bounds, workl, workl(np+1))
         if (nev .eq. nev0+1) numcnv = nev0+1
!
!        %-------------------%
!        | Convergence test. |
!        %-------------------%
!
         call dcopy  (nev, bounds(np+1), 1, workl(2*np+1), 1)
         call dnconv  (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1),
     &        tol, nconv)
!
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call pivout (comm, logfil, 4, kp, ndigit,
     &                  '_naup2: NEV, NP, NUMCNV, NCONV are')
            call pdvout  (comm, logfil, kplusp, ritzr, ndigit,
     &           '_naup2: Real part of the eigenvalues of H')
            call pdvout  (comm, logfil, kplusp, ritzi, ndigit,
     &           '_naup2: Imaginary part of the eigenvalues of H')
            call pdvout  (comm, logfil, kplusp, bounds, ndigit,
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
!
!        %---------------------------------------------------------%
!        | Count the number of unwanted Ritz values that have zero |
!        | Ritz estimates. If any Ritz estimates are equal to zero |
!        | then a leading block of H of order equal to at least    |
!        | the number of Ritz values with zero Ritz estimates has  |
!        | split off. None of these Ritz values may be removed by  |
!        | shifting. Decrease NP the number of shifts to apply. If |
!        | no shifts may be applied, then prepare to exit          |
!        %---------------------------------------------------------%
!
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
!
         if ( (nconv .ge. numcnv) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
!
            if (msglvl .gt. 4) then
               call dvout (logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Real part of the eig computed by _neigh:')
               call dvout (logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Imag part of the eig computed by _neigh:')
               call dvout (logfil, kplusp, workl(kplusp**2+kplusp*2+1),
     &                     ndigit,
     &             '_naup2: Ritz estimates computed by _neigh:')
            end if
!
!           %------------------------------------------------%
!           | Prepare to exit. Put the converged Ritz values |
!           | and corresponding bounds in RITZ(1:NCONV) and  |
!           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
!           | careful when NCONV > NP                        |
!           %------------------------------------------------%
!
!           %------------------------------------------%
!           |  Use h( 3,1 ) as storage to communicate  |
!           |  rnorm to _neupd if needed               |
!           %------------------------------------------%

            h(3,1) = rnorm
!
!           %----------------------------------------------%
!           | To be consistent with dngets , we first do a  |
!           | pre-processing sort in order to keep complex |
!           | conjugate pairs together.  This is similar   |
!           | to the pre-processing sort used in dngets     |
!           | except that the sort is done in the opposite |
!           | order.                                       |
!           %----------------------------------------------%
!
            if (which .eq. 'LM') wprime = 'SR'
            if (which .eq. 'SM') wprime = 'LR'
            if (which .eq. 'LR') wprime = 'SM'
            if (which .eq. 'SR') wprime = 'LM'
            if (which .eq. 'LI') wprime = 'SM'
            if (which .eq. 'SI') wprime = 'LM'
!
            call dsortc  (wprime, .true., kplusp, ritzr, ritzi, bounds)
!
!           %----------------------------------------------%
!           | Now sort Ritz values so that converged Ritz  |
!           | values appear within the first NEV locations |
!           | of ritzr, ritzi and bounds, and the most     |
!           | desired one appears at the front.            |
!           %----------------------------------------------%
!
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
!
            call dsortc (wprime, .true., kplusp, ritzr, ritzi, bounds)
!
!           %--------------------------------------------------%
!           | Scale the Ritz estimate of each Ritz value       |
!           | by 1 / max(eps23,magnitude of the Ritz value).   |
!           %--------------------------------------------------%
!
            do 35 j = 1, numcnv
                temp = max(eps23,dlapy2 (ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
!
!           %----------------------------------------------------%
!           | Sort the Ritz values according to the scaled Ritz  |
!           | esitmates.  This will push all the converged ones  |
!           | towards the front of ritzr, ritzi, bounds          |
!           | (in the case when NCONV < NEV.)                    |
!           %----------------------------------------------------%
!
            wprime = 'LR'
            call dsortc (wprime, .true., numcnv, bounds, ritzr, ritzi)
!
!           %----------------------------------------------%
!           | Scale the Ritz estimate back to its original |
!           | value.                                       |
!           %----------------------------------------------%
!
            do 40 j = 1, numcnv
                temp = max(eps23, dlapy2 (ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)*temp
 40         continue
!
!           %------------------------------------------------%
!           | Sort the converged Ritz values again so that   |
!           | the "threshold" value appears at the front of  |
!           | ritzr, ritzi and bound.                        |
!           %------------------------------------------------%
!
            call dsortc (which, .true., nconv, ritzr, ritzi, bounds)
!
!
            if (msglvl .gt. 1) then
               call dvout  (logfil, kplusp, ritzr, ndigit,
     &            '_naup2: Sorted real part of the eigenvalues')
               call dvout  (logfil, kplusp, ritzi, ndigit,
     &            '_naup2: Sorted imaginary part of the eigenvalues')
               call dvout  (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if
!
!           %------------------------------------%
!           | Max iterations have been exceeded. |
!           %------------------------------------%
!
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
!
!           %---------------------%
!           | No shifts to apply. |
!           %---------------------%
!
            if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
!
            np = nconv
            go to 1100
!
         else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
!
!           %-------------------------------------------------%
!           | Do not have all the requested eigenvalues yet.  |
!           | To prevent possible stagnation, adjust the size |
!           | of NEV.                                         |
!           %-------------------------------------------------%
!
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
!
!           %---------------------------------------%
!           | If the size of NEV was just increased |
!           | resort the eigenvalues.               |
!           %---------------------------------------%
!
            if (nevbef .lt. nev)
     &         call pdngets (comm, ishift, which, nev, np, ritzr, ritzi,
     &                      bounds, workl, workl(np+1))
!
         end if
!
         if (msglvl .gt. 0) then
            call pivout (comm, logfil, 1, [nconv], ndigit,
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call pivout (comm, logfil, 2, kp, ndigit,
     &              '_naup2: NEV and NP are')
               call pdvout  (comm, logfil, nev, ritzr(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- real part')
               call pdvout  (comm, logfil, nev, ritzi(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- imag part')
               call pdvout  (comm, logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
!
         if (ishift .eq. 0) then
!
!           %-------------------------------------------------------%
!           | User specified shifts: reverse communication to       |
!           | compute the shifts. They are returned in the first    |
!           | 2*NP locations of WORKL.                              |
!           %-------------------------------------------------------%
!
            ushift = .true.
            ido = 3
            go to 9000
         end if
!
   50    continue
!
!        %------------------------------------%
!        | Back from reverse communication;   |
!        | User specified shifts are returned |
!        | in WORKL(1:2*NP)                   |
!        %------------------------------------%
!
         ushift = .false.
!
         if ( ishift .eq. 0 ) then
!
!            %----------------------------------%
!            | Move the NP shifts from WORKL to |
!            | RITZR, RITZI to free up WORKL    |
!            | for non-exact shift case.        |
!            %----------------------------------%
!
             call dcopy  (np, workl,       1, ritzr, 1)
             call dcopy  (np, workl(np+1), 1, ritzi, 1)
         end if
!
         if (msglvl .gt. 2) then
            call pivout (comm, logfil, 1, [np], ndigit,
     &                  '_naup2: The number of shifts to apply ')
            call pdvout  (comm, logfil, np, ritzr, ndigit,
     &                  '_naup2: Real part of the shifts')
            call pdvout  (comm, logfil, np, ritzi, ndigit,
     &                  '_naup2: Imaginary part of the shifts')
            if ( ishift .eq. 1 )
     &          call pdvout  (comm, logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if
!
!        %---------------------------------------------------------%
!        | Apply the NP implicit shifts by QR bulge chasing.       |
!        | Each shift is applied to the whole upper Hessenberg     |
!        | matrix H.                                               |
!        | The first 2*N locations of WORKD are used as workspace. |
!        %---------------------------------------------------------%
!
         call pdnapps  (comm, n, nev, np, ritzr, ritzi, v, ldv,
     &                 h, ldh, resid, q, ldq, workl, workd)
!
!        %---------------------------------------------%
!        | Compute the B-norm of the updated residual. |
!        | Keep B*RESID in WORKD(1:N) to be used in    |
!        | the first step of the next call to pdnaitr . |
!        %---------------------------------------------%
!
         cnorm = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy  (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*RESID |
!           %----------------------------------%
!
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy  (n, resid, 1, workd, 1)
         end if
!
  100    continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(1:N) := B*RESID            |
!        %----------------------------------%
!
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         endif
!
         if (bmat .eq. 'G') then
            rnorm_buf = ddot  (n, resid, 1, workd, 1)
            call MPI_ALLREDUCE( [rnorm_buf], buf2, 1,
     &                MPI_DOUBLE_PRECISION , MPI_SUM, comm, ierr )
            rnorm = sqrt(abs(buf2(1)))
         else if (bmat .eq. 'I') then
            rnorm = pdnorm2 ( comm, n, resid, 1 )
         end if
         cnorm = .false.
!
         if (msglvl .gt. 2) then
            call pdvout  (comm, logfil, 1, [rnorm], ndigit,
     &      '_naup2: B-norm of residual for compressed factorization')
            call pdmout  (comm, logfil, nev, nev, h, ldh, ndigit,
     &        '_naup2: Compressed upper Hessenberg matrix H')
         end if
!
      go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 1100 continue
!
      mxiter = iter
      nev = numcnv
!
 1200 continue
      ido = 99
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      call arscnd (t1)
      tnaup2 = t1 - t0
!
 9000 continue
!
!     %----------------%
!     | End of pdnaup2  |
!     %----------------%
!
      return
      end
