!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: pssaup2
!
! Message Passing Layer: MPI
!
!\Description:
!  Intermediate level interface called by pssaupd.
!
!\Usage:
!  call pssaup2
!     ( COMM, IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL,
!       IPNTR, WORKD, INFO )
!
!\Arguments
!
!  COMM, IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in pssaupd.
!  MODE, ISHIFT, MXITER: see the definition of IPARAM in pssaupd.
!
!  NP      Integer.  (INPUT/OUTPUT)
!          Contains the number of implicit shifts to apply during
!          each Arnoldi/Lanczos iteration.
!          If ISHIFT=1, NP is adjusted dynamically at each iteration
!          to accelerate convergence and prevent stagnation.
!          This is also roughly equal to the number of matrix-vector
!          products (involving the operator OP) per Arnoldi iteration.
!          The logic for adjusting is contained within the current
!          subroutine.
!          If ISHIFT=0, NP is the number of shifts the user needs
!          to provide via reverse communication. 0 < NP < NCV-NEV.
!          NP may be less than NCV-NEV since a leading block of the current
!          upper Tridiagonal matrix has split off and contains "unwanted"
!          Ritz values.
!          Upon termination of the IRA iteration, NP contains the number
!          of "converged" wanted Ritz values.
!
!  IUPD    Integer.  (INPUT)
!          IUPD .EQ. 0: use explicit restart instead implicit update.
!          IUPD .NE. 0: use implicit update.
!
!  V       Real N by (NEV+NP) array.  (INPUT/OUTPUT)
!          The Lanczos basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (NEV+NP) by 2 array.  (OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          The subdiagonal is stored in the first column of H starting
!          at H(2,1).  The main diagonal is stored in the second column
!          of H starting at H(1,2). If pssaup2 converges store the
!          B-norm of the final residual vector in H(1,1).
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZ    Real array of length NEV+NP.  (OUTPUT)
!          RITZ(1:NEV) contains the computed Ritz values of OP.
!
!  BOUNDS  Real array of length NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
!
!  Q       Real (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
!          Private (replicated) work array used to accumulate the
!          rotation in the shift application step.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Real array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  It is used in the computation of the
!          tridiagonal eigenvalue problem, the calculation and
!          application of the shifts and convergence checking.
!          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
!          of WORKL are used in reverse communication to hold the user
!          supplied shifts.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD for
!          vectors used by the Lanczos iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in one of
!                    the spectral transformation modes.  X is the current
!                    operand.
!          -------------------------------------------------------------
!
!  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Lanczos iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration !!!!!!!!!!
!          See Data Distribution Note in pssaupd.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =     0: Normal return.
!          =     1: All possible eigenvalues of OP has been found.
!                   NP returns the size of the invariant subspace
!                   spanning the operator OP.
!          =     2: No shifts could be applied.
!          =    -8: Error return from trid. eigenvalue calculation;
!                   This should never happen.
!          =    -9: Starting vector is zero.
!          = -9999: Could not build an Lanczos factorization.
!                   Size that was built in returned in NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Routines called:
!     psgetv0  Parallel ARPACK initial vector generation routine.
!     pssaitr  Parallel ARPACK Lanczos factorization routine.
!     pssapps  Parallel ARPACK application of implicit shifts routine.
!     ssconv   ARPACK convergence of Ritz values routine.
!     psseigt  Parallel ARPACK compute Ritz values and error bounds routine.
!     pssgets  Parallel ARPACK reorder Ritz values and error bounds routine.
!     ssortr   ARPACK sorting routine.
!     sstrqb   ARPACK routine that computes all eigenvalues and the
!              last component of the eigenvectors of a symmetric
!              tridiagonal matrix using the implicit QL or QR method.
!     pivout   Parallel ARPACK utility routine that prints integers.
!     arscnd   ARPACK utility routine for timing.
!     psvout   Parallel ARPACK utility routine that prints vectors.
!     pslamch10  ScaLAPACK routine that determines machine constants.
!     scopy    Level 1 BLAS that copies one vector to another.
!     sdot     Level 1 BLAS that computes the scalar product of two vectors.
!     psnorm2  Parallel version of Level 1 BLAS that computes the norm of a vector.
!     sscal    Level 1 BLAS that scales a vector.
!     sswap    Level 1 BLAS that swaps two vectors.
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
!     Starting Point: Serial Code FILE: saup2.F   SID: 2.4
!
!\SCCS Information:
! FILE: saup2.F   SID: 1.5   DATE OF SID: 05/20/98
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pssaup2
     &   ( comm, ido, bmat, n, which, nev, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, h, ldh, ritz, bounds,
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
      integer    ido, info, ishift, iupd, ldh, ldq, ldv, mxiter,
     &           n, mode, nev, np
      Real
     &           tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Real
     &           bounds(nev+np), h(ldh,2), q(ldq,nev+np), resid(n),
     &           ritz(nev+np), v(ldv,nev+np), workd(3*n),
     &           workl(3*(nev+np))
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift
      integer    ierr, iter, j, kplusp, msglvl, nconv, nevbef, nev0,
     &           np0, nptemp, nevd2, nevm2, kp(3)
      Real
     &           rnorm, temp, eps23, buf2(1)
      save       cnorm, getv0, initv, update, ushift,
     &           iter, kplusp, msglvl, nconv, nev0, np0,
     &           rnorm, eps23
!
      Real
     &           rnorm_buf
!
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy, psgetv0, pssaitr, sscal, ssconv,
     &           psseigt, pssgets, pssapps,
     &           ssortr, psvout, pivout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           sdot, psnorm2, pslamch10
      external   sdot, psnorm2, pslamch10
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min, abs, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido .eq. 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call arscnd (t0)
         msglvl = msaup2
!
!        %---------------------------------%
!        | Set machine dependent constant. |
!        %---------------------------------%
!
         eps23 = pslamch10(comm, 'Epsilon-Machine')
         eps23 = eps23**(2.0/3.0)
!
!        %-------------------------------------%
!        | nev0 and np0 are integer variables  |
!        | hold the initial values of NEV & NP |
!        %-------------------------------------%
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
         kplusp = nev0 + np0
         nconv  = 0
         iter   = 0
!
!        %---------------------------------------------%
!        | Set flags for computing the first NEV steps |
!        | of the Lanczos factorization.               |
!        %---------------------------------------------%
!
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
!
         if (info .ne. 0) then
!
!        %--------------------------------------------%
!        | User provides the initial residual vector. |
!        %--------------------------------------------%
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
         call psgetv0 ( comm, ido, bmat, 1, initv, n, 1, v, ldv,
     &                  resid, rnorm, ipntr, workd, workl, info)
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
            go to 1200
         end if
         getv0 = .false.
         ido  = 0
      end if
!
!     %------------------------------------------------------------%
!     | Back from reverse communication: continue with update step |
!     %------------------------------------------------------------%
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
!     | Compute the first NEV steps of the Lanczos factorization |
!     %----------------------------------------------------------%
!
      call pssaitr (comm, ido, bmat, n, 0, nev0, mode,
     &              resid, rnorm, v, ldv, h, ldh, ipntr,
     &              workd, workl, info)
!
!     %---------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication  |
!     | to compute operations involving OP and possibly B |
!     %---------------------------------------------------%
!
      if (ido .ne. 99) go to 9000
!
      if (info .gt. 0) then
!
!        %-----------------------------------------------------%
!        | pssaitr was unable to build an Lanczos factorization|
!        | of length NEV0. INFO is returned with the size of   |
!        | the factorization built. Exit main loop.            |
!        %-----------------------------------------------------%
!
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |
!     |           Each iteration implicitly restarts the Lanczos     |
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
     &           '_saup2: **** Start of major iteration number ****')
         end if
         if (msglvl .gt. 1) then
            call pivout (comm, logfil, 1, [nev], ndigit,
     &     '_saup2: The length of the current Lanczos factorization')
            call pivout (comm, logfil, 1, [np], ndigit,
     &           '_saup2: Extend the Lanczos factorization by')
         end if
!
!        %------------------------------------------------------------%
!        | Compute NP additional steps of the Lanczos factorization.  |
!        %------------------------------------------------------------%
!
         ido = 0
   20    continue
         update = .true.
!
         call pssaitr (comm, ido, bmat, n, nev, np, mode,
     &                 resid, rnorm, v, ldv, h, ldh, ipntr,
     &                 workd, workl, info)
!
!        %---------------------------------------------------%
!        | ido .ne. 99 implies use of reverse communication  |
!        | to compute operations involving OP and possibly B |
!        %---------------------------------------------------%
!
         if (ido .ne. 99) go to 9000
!
         if (info .gt. 0) then
!
!           %-----------------------------------------------------%
!           | pssaitr was unable to build an Lanczos factorization|
!           | of length NEV0+NP0. INFO is returned with the size  |
!           | of the factorization built. Exit main loop.         |
!           %-----------------------------------------------------%
!
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
!
         if (msglvl .gt. 1) then
            call psvout (comm, logfil, 1, [rnorm], ndigit,
     &           '_saup2: Current B-norm of residual for factorization')
         end if
!
!        %--------------------------------------------------------%
!        | Compute the eigenvalues and corresponding error bounds |
!        | of the current symmetric tridiagonal matrix.           |
!        %--------------------------------------------------------%
!
         call psseigt ( comm, rnorm, kplusp, h, ldh, ritz, bounds,
     &                  workl, ierr)
!
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
!
!        %----------------------------------------------------%
!        | Make a copy of eigenvalues and corresponding error |
!        | bounds obtained from _seigt.                       |
!        %----------------------------------------------------%
!
         call scopy(kplusp, ritz, 1, workl(kplusp+1), 1)
         call scopy(kplusp, bounds, 1, workl(2*kplusp+1), 1)
!
!        %---------------------------------------------------%
!        | Select the wanted Ritz values and their bounds    |
!        | to be used in the convergence test.               |
!        | The selection is based on the requested number of |
!        | eigenvalues instead of the current NEV and NP to  |
!        | prevent possible misconvergence.                  |
!        | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |
!        | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |
!        %---------------------------------------------------%
!
         nev = nev0
         np = np0
         call pssgets ( comm, ishift, which, nev, np, ritz,
     &                  bounds, workl)
!
!        %-------------------%
!        | Convergence test. |
!        %-------------------%
!
         call scopy (nev, bounds(np+1), 1, workl(np+1), 1)
         call ssconv (nev, ritz(np+1), workl(np+1), tol, nconv)
!
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call pivout (comm, logfil, 3, kp, ndigit,
     &                  '_saup2: NEV, NP, NCONV are')
            call psvout (comm, logfil, kplusp, ritz, ndigit,
     &           '_saup2: The eigenvalues of H')
            call psvout (comm, logfil, kplusp, bounds, ndigit,
     &          '_saup2: Ritz estimates of the current NCV Ritz values')
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
         if ( (nconv .ge. nev0) .or.
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
!
!           %------------------------------------------------%
!           | Prepare to exit. Put the converged Ritz values |
!           | and corresponding bounds in RITZ(1:NCONV) and  |
!           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
!           | careful when NCONV > NP since we don't want to |
!           | swap overlapping locations.                    |
!           %------------------------------------------------%
!
            if (which .eq. 'BE') then
!
!              %-----------------------------------------------------%
!              | Both ends of the spectrum are requested.            |
!              | Sort the eigenvalues into algebraically decreasing  |
!              | order first then swap low end of the spectrum next  |
!              | to high end in appropriate locations.               |
!              | NOTE: when np < floor(nev/2) be careful not to swap |
!              | overlapping locations.                              |
!              %-----------------------------------------------------%
!
               wprime = 'SA'
               call ssortr (wprime, .true., kplusp, ritz, bounds)
               nevd2 = nev0 / 2
               nevm2 = nev0 - nevd2
               if ( nev .gt. 1 ) then
                  call sswap ( min(nevd2,np), ritz(nevm2+1), 1,
     &                 ritz( max(kplusp-nevd2+1,kplusp-np+1) ), 1)
                  call sswap ( min(nevd2,np), bounds(nevm2+1), 1,
     &                 bounds( max(kplusp-nevd2+1,kplusp-np+1)), 1)
               end if
!
            else
!
!              %--------------------------------------------------%
!              | LM, SM, LA, SA case.                             |
!              | Sort the eigenvalues of H into the an order that |
!              | is opposite to WHICH, and apply the resulting    |
!              | order to BOUNDS.  The eigenvalues are sorted so  |
!              | that the wanted part are always within the first |
!              | NEV locations.                                   |
!              %--------------------------------------------------%
!
               if (which .eq. 'LM') wprime = 'SM'
               if (which .eq. 'SM') wprime = 'LM'
               if (which .eq. 'LA') wprime = 'SA'
               if (which .eq. 'SA') wprime = 'LA'
!
               call ssortr (wprime, .true., kplusp, ritz, bounds)
!
            end if
!
!           %--------------------------------------------------%
!           | Scale the Ritz estimate of each Ritz value       |
!           | by 1 / max(eps23,magnitude of the Ritz value).   |
!           %--------------------------------------------------%
!
            do 35 j = 1, nev0
               temp = max( eps23, abs(ritz(j)) )
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
            wprime = 'LA'
            call ssortr(wprime, .true., nev0, bounds, ritz)
!
!           %----------------------------------------------%
!           | Scale the Ritz estimate back to its original |
!           | value.                                       |
!           %----------------------------------------------%
!
            do 40 j = 1, nev0
                temp = max( eps23, abs(ritz(j)) )
                bounds(j) = bounds(j)*temp
 40         continue
!
!           %--------------------------------------------------%
!           | Sort the "converged" Ritz values again so that   |
!           | the "threshold" values and their associated Ritz |
!           | estimates appear at the appropriate position in  |
!           | ritz and bound.                                  |
!           %--------------------------------------------------%
!
            if (which .eq. 'BE') then
!
!              %------------------------------------------------%
!              | Sort the "converged" Ritz values in increasing |
!              | order.  The "threshold" values are in the      |
!              | middle.                                        |
!              %------------------------------------------------%
!
               wprime = 'LA'
               call ssortr(wprime, .true., nconv, ritz, bounds)
!
            else
!
!              %----------------------------------------------%
!              | In LM, SM, LA, SA case, sort the "converged" |
!              | Ritz values according to WHICH so that the   |
!              | "threshold" value appears at the front of    |
!              | ritz.                                        |
!              %----------------------------------------------%

               call ssortr(which, .true., nconv, ritz, bounds)
!
            end if
!
!           %------------------------------------------%
!           |  Use h( 1,1 ) as storage to communicate  |
!           |  rnorm to _seupd if needed               |
!           %------------------------------------------%
!
            h(1,1) = rnorm
!
            if (msglvl .gt. 1) then
               call psvout (comm, logfil, kplusp, ritz, ndigit,
     &            '_saup2: Sorted Ritz values.')
               call psvout (comm, logfil, kplusp, bounds, ndigit,
     &            '_saup2: Sorted ritz estimates.')
            end if
!
!           %------------------------------------%
!           | Max iterations have been exceeded. |
!           %------------------------------------%
!
            if (iter .gt. mxiter .and. nconv .lt. nev) info = 1
!
!           %---------------------%
!           | No shifts to apply. |
!           %---------------------%
!
            if (np .eq. 0 .and. nconv .lt. nev0) info = 2
!
            np = nconv
            go to 1100
!
         else if (nconv .lt. nev .and. ishift .eq. 1) then
!
!           %---------------------------------------------------%
!           | Do not have all the requested eigenvalues yet.    |
!           | To prevent possible stagnation, adjust the number |
!           | of Ritz values and the shifts.                    |
!           %---------------------------------------------------%
!
            nevbef = nev
            nev = nev + min (nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 2) then
               nev = 2
            end if
            np  = kplusp - nev
!
!           %---------------------------------------%
!           | If the size of NEV was just increased |
!           | resort the eigenvalues.               |
!           %---------------------------------------%
!
            if (nevbef .lt. nev)
     &         call pssgets ( comm, ishift, which, nev, np,
     &                        ritz, bounds, workl)
!
         end if
!
         if (msglvl .gt. 0) then
            call pivout (comm, logfil, 1, [nconv], ndigit,
     &           '_saup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call pivout (comm, logfil, 2, kp, ndigit,
     &              '_saup2: NEV and NP .')
               call psvout (comm, logfil, nev, ritz(np+1), ndigit,
     &              '_saup2: "wanted" Ritz values.')
               call psvout (comm, logfil, nev, bounds(np+1), ndigit,
     &              '_saup2: Ritz estimates of the "wanted" values ')
            end if
         end if
!
         if (ishift .eq. 0) then
!
!           %-----------------------------------------------------%
!           | User specified shifts: reverse communication to     |
!           | compute the shifts. They are returned in the first  |
!           | NP locations of WORKL.                              |
!           %-----------------------------------------------------%
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
!        | in WORKL(1:NP)                     |
!        %------------------------------------%
!
         ushift = .false.
!
!
!        %---------------------------------------------------------%
!        | Move the NP shifts to the first NP locations of RITZ to |
!        | free up WORKL.  This is for the non-exact shift case;   |
!        | in the exact shift case, pssgets already handles this.  |
!        %---------------------------------------------------------%
!
         if (ishift .eq. 0) call scopy (np, workl, 1, ritz, 1)
!
         if (msglvl .gt. 2) then
            call pivout (comm, logfil, 1, [np], ndigit,
     &                  '_saup2: The number of shifts to apply ')
            call psvout (comm, logfil, np, workl, ndigit,
     &                  '_saup2: shifts selected')
            if (ishift .eq. 1) then
               call psvout (comm, logfil, np, bounds, ndigit,
     &                  '_saup2: corresponding Ritz estimates')
             end if
         end if
!
!        %---------------------------------------------------------%
!        | Apply the NP0 implicit shifts by QR bulge chasing.      |
!        | Each shift is applied to the entire tridiagonal matrix. |
!        | The first 2*N locations of WORKD are used as workspace. |
!        | After pssapps is done, we have a Lanczos                |
!        | factorization of length NEV.                            |
!        %---------------------------------------------------------%
!
         call pssapps ( comm, n, nev, np, ritz, v, ldv, h, ldh, resid,
     &                  q, ldq, workd)
!
!        %---------------------------------------------%
!        | Compute the B-norm of the updated residual. |
!        | Keep B*RESID in WORKD(1:N) to be used in    |
!        | the first step of the next call to pssaitr. |
!        %---------------------------------------------%
!
         cnorm = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(n+1), 1)
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
            call scopy (n, resid, 1, workd, 1)
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
         end if
!
         if (bmat .eq. 'G') then
            rnorm_buf = sdot (n, resid, 1, workd, 1)
            call MPI_ALLREDUCE( [rnorm_buf], buf2, 1,
     &                MPI_REAL, MPI_SUM, comm, ierr )
            rnorm = sqrt(abs(buf2(1)))
         else if (bmat .eq. 'I') then
            rnorm = psnorm2( comm, n, resid, 1 )
         end if
         cnorm = .false.
  130    continue
!
         if (msglvl .gt. 2) then
            call psvout (comm, logfil, 1, [rnorm], ndigit,
     &      '_saup2: B-norm of residual for NEV factorization')
            call psvout (comm, logfil, nev, h(1,2), ndigit,
     &           '_saup2: main diagonal of compressed H matrix')
            call psvout (comm, logfil, nev-1, h(2,1), ndigit,
     &           '_saup2: subdiagonal of compressed H matrix')
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
      nev = nconv
!
 1200 continue
      ido = 99
!
!     %------------%
!     | Error exit |
!     %------------%
!
      call arscnd (t1)
      tsaup2 = t1 - t0
!
 9000 continue
      return
!
!     %----------------%
!     | End of pssaup2 |
!     %----------------%
!
      end
