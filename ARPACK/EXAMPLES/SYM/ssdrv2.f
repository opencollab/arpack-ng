      program ssdrv2 
c
c     Program to illustrate the idea of reverse communication
c     in shift and invert mode for a standard symmetric eigenvalue
c     problem.  The following program uses the two LAPACK subroutines 
c     sgttrf.f and sgttrs.f to factor and solve a tridiagonal system of 
c     equations.
c
c     We implement example two of ex-sym.doc in DOCUMENTS directory
c
c\Example-2
c     ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
c         where A is derived from the central difference discretization
c         of the 1-dimensional Laplacian on [0,1]  with zero Dirichlet
c         boundary condition.
c     ... OP = (inv[A - sigma*I]) and  B = I.
c     ... Use mode 3 of SSAUPD.
c
c\BeginLib
c
c\Routines called:
c     ssaupd  ARPACK reverse communication interface routine.
c     sseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     sgttrf  LAPACK tridiagonal factorization routine.
c     sgttrs  LAPACK tridiagonal solve routine.
c     saxpy   saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     snrm2   Level 1 BLAS that computes the norm of a vector.     
c     av      Matrix vector multiplication routine that computes A*x.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: sdrv2.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c----------------------------------------------------------------------
c
c     %-----------------------------%
c     | Define leading dimensions   |
c     | for all arrays.             |
c     | MAXN:   Maximum dimension   |
c     |         of the A allowed.   |
c     | MAXNEV: Maximum NEV allowed |
c     | MAXNCV: Maximum NCV allowed |
c     %-----------------------------%
c
      integer          maxn, maxnev, maxncv, ldv
      parameter        (maxn=256, maxnev=10, maxncv=25, 
     &                 ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Real
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxn), d(maxncv,2), resid(maxn), 
     &                 ad(maxn), adl(maxn), adu(maxn), adu2(maxn),
     &                 ax(maxn)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11), ipiv(maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, j, ierr,
     &                 nconv, maxitr, ishfts, mode
      logical          rvec
      Real
     &                 sigma, tol, h2
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real           
     &                 zero, one, two 
      parameter        (zero = 0.0E+0, one = 1.0E+0, 
     &                  two = 2.0E+0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real
     &                 snrm2
      external         saxpy, snrm2, sgttrf, sgttrs
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic        abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %----------------------------------------------------%
c     | The number N is the dimension of the matrix.  A    |
c     | standard eigenvalue problem is solved (BMAT = 'I'. |
c     | NEV is the number of eigenvalues (closest to       |
c     | SIGMA) to be approximated.  Since the shift-invert |
c     | mode is used, WHICH is set to 'LM'.  The user can  |
c     | modify NEV, NCV, SIGMA to solve problems of        | 
c     | different sizes, and to get different parts of the |
c     | spectrum.  However, The following conditions must  |
c     | be satisfied:                                      |
c     |                   N <= MAXN,                       | 
c     |                 NEV <= MAXNEV,                     |
c     |             NEV + 1 <= NCV <= MAXNCV               | 
c     %----------------------------------------------------% 
c
      n = 100
      nev = 4
      ncv = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SDRV2: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV2: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV2: NCV is greater than MAXNCV '
         go to 9000
      end if
c
      bmat = 'I'
      which = 'LM'
      sigma = zero 
c
c     %--------------------------------------------------%
c     | The work array WORKL is used in SSAUPD as        |
c     | workspace.  Its dimension LWORKL is set as       |
c     | illustrated below.  The parameter TOL determines |
c     | the stopping criterion.  If TOL<=0, machine      |
c     | precision is used.  The variable IDO is used for |
c     | reverse communication and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is |
c     | generated in SSAUPD to start the Arnoldi         |
c     | iteration.                                       |
c     %--------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol = zero 
      ido = 0
      info = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 3 of SSAUPD is used     |
c     | (IPARAM(7) = 3).  All these options may be        |
c     | changed by the user. For details, see the         |
c     | documentation in SSAUPD.                          |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 3
c
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
c
c     %-----------------------------------------------------%
c     | Call LAPACK routine to factor (A-SIGMA*I), where A  |
c     | is the 1-d Laplacian.                               | 
c     %-----------------------------------------------------%
c
      h2 = one / real((n+1)*(n+1))
      do 20 j=1,n
         ad(j) = two / h2 - sigma
         adl(j) = -one / h2
 20   continue 
      call scopy (n, adl, 1, adu, 1)
      call sgttrf (n, adl, ad, adu, adu2, ipiv, ierr)
      if (ierr .ne. 0) then 
         print *, ' ' 
         print *, ' Error with _gttrf in SDRV2.'
         print *, ' '
         go to 9000
      end if
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine SSAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call ssaupd ( ido, bmat, n, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %----------------------------------------%
c           | Perform y <-- OP*x = inv[A-sigma*I]*x. |
c           | The user only need the linear system   |
c           | solver here that takes workd(ipntr(1)) |
c           | as input, and returns the result to    |
c           | workd(ipntr(2)).                       |
c           %----------------------------------------%
c
            call scopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
c
            call sgttrs ('Notranspose', n, 1, adl, ad, adu, adu2, ipiv,
     &                   workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print *, ' '
               print *, ' Error with _gttrs in _SDRV2. '
               print *, ' ' 
               go to 9000
            end if
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %----------------------------%
c        | Error message.  Check the  |
c        | documentation in SSAUPD    |
c        %----------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ',info
         print *, ' Check documentation of _saupd '
         print *, ' ' 
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using SSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call sseupd ( rvec, 'All', select, d, v, ldv, sigma,
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &        iparam, ipntr, workd, workl, lworkl, ierr )
c
c        %----------------------------------------------%
c        | Eigenvalues are returned in the first column |
c        | of the two dimensional array D and the       |
c        | corresponding eigenvectors are returned in   |
c        | the first NEV columns of the two dimensional |
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%

         if ( ierr .ne. 0 ) then 
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of SSEUPD. |
c           %------------------------------------%
c
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd '
            print *, ' '
c
         else 
c
            nconv =  iparam(5)
            do 30 j=1, nconv
c
c              %---------------------------%
c              | Compute the residual norm |
c              |                           |
c              |   ||  A*x - lambda*x ||   |
c              |                           |
c              | for the NCONV accurately  |
c              | computed eigenvalues and  |
c              | eigenvectors.  (iparam(5) |
c              | indicates how many are    |
c              | accurate to the requested |
c              | tolerance)                |
c              %---------------------------%
c
               call av(n, v(1,j), ax)
               call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
               d(j,2) = snrm2(n, ax, 1)
               d(j,2) = d(j,2) / abs(d(j,1))
c
 30         continue
c
c           %-------------------------------%
c           | Display computed residuals    |
c           %-------------------------------%
c
            call smout(6, nconv, 2, d, maxncv, -6,
     &           'Ritz values and relative residuals')
         end if
c
c        %------------------------------------------%
c        | Print additional convergence information |
c        %------------------------------------------%
c
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit',
     &               ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
c
         print *, ' '
         print *, ' _SDRV2 '
         print *, ' ====== '
         print *, ' ' 
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',
     &              nconv 
         print *, ' The number of Implicit Arnoldi update',
     &            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
c
      end if
c
c     %---------------------------%
c     | Done with program ssdrv2. |
c     %---------------------------%
c
 9000 continue
c
      end
c
c------------------------------------------------------------------------
c     Matrix vector subroutine
c     where the matrix is the 1 dimensional discrete Laplacian on
c     the interval [0,1] with zero Dirichlet boundary condition.
c
      subroutine av (n, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), one, two, h2
      parameter         (one = 1.0E+0,  two = 2.0E+0)

c
      w(1) =  two*v(1) - v(2)
      do 100 j = 2,n-1
         w(j) = - v(j-1) + two*v(j) - v(j+1) 
  100 continue
      j = n
      w(j) = - v(j-1) + two*v(j) 
c
c     Scale the vector w by (1 / h^2).
c
      h2 = one / real((n+1)*(n+1))
      call sscal(n, one/h2, w, 1)
      return
      end
