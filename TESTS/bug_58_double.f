      program dndrv2
      implicit none
c
c     Test program to show that an eigenvector for a modified, singular
c     9x9 identity matrix is a column of NaNs. See issue #58. The problem
c     is in the purification stage of dneupd.
c
c     The shift sigma is the real number -1.0d0.
c
c     OP = inv[A-sigma*I] and  B = I.
c
c     Use mode 3 of DNAUPD.
c
c\BeginLib
c
c\Routines called:
c     dnaupd  ARPACK reverse communication interface routine.
c     dneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dgttrf  LAPACK tridiagonal factorization routine.
c     dgttrs  LAPACK tridiagonal solve routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
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
c FILE: ndrv2.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c--------------------------------------------------------------------------
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
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=9, maxnev=9, maxncv=9,
     &                   ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14), ipiv(maxn)
      logical           select(maxncv)
      Double precision
     &                  ax(maxn), d(maxncv,3), resid(maxn),
     &                  v(ldv, maxncv), workd(3*maxn),
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv),
     &                  dd(maxn), dl(maxn), du(maxn),
     &                  du2(maxn), a(maxn,maxn), c(maxn,maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j, i, k,
     &                  ierr, nconv, maxitr, ishfts, mode, nnz
      Double precision
     &                  tol, h, s,
     &                  sigmar, sigmai, s1, s2, s3
      logical           first, rvec
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                   one, zero, two, rho
      common            /convct/ rho
      parameter         (one = 1.0D+0, zero = 0.0D+0,
     &                   two = 2.0D+0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision
     &                  ddot, dnrm2, dlapy2
      external          dgttrf, dgttrs, ddot, dnrm2, dlapy2
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic         abs
c
c     %-----------------------%
c     | Executable statements |
c     %-----------------------%
c
c     %--------------------------------------------------%
c     | The number N is the dimension of the matrix.  A  |
c     | standard eigenvalue problem is solved (BMAT =    |
c     | 'I').  NEV is the number of eigenvalues (closest |
c     | to the shift SIGMAR) to be approximated.  Since  |
c     | the shift-invert mode is used, WHICH is set to   |
c     | 'LM'.  The user can modify NEV, NCV, SIGMAR to   |
c     | solve problems of different sizes, and to get    |
c     | different parts of the spectrum.  However, The   |
c     | following conditions must be satisfied:          |
c     |                 N <= MAXN,                       |
c     |               NEV <= MAXNEV,                     |
c     |           NEV + 2 <= NCV <= MAXNCV               |
c     %--------------------------------------------------%
c
      nev   = 4
      ncv   = 8
      do i = 1,maxn
         do j = 1,maxn
            a(i,j) = 0.0d0
            c(i,j) = 0.0d0
         end do
      end do
      n = 9
      do i = 1,n
         a(i,i) = 1.0d0
         c(i,i) = 1.0d0
      end do
      a(1,1) = 0.0d0
      c(1,1) = 0.0d0
      a(1,n) = 1.0d0
      c(1,n) = 1.0d0
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV2: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV2: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV2: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'LM'
      sigmar = -1.0D+0
      sigmai = 0.0D+0
c
c     %----------------------------------------------------%
c     | Construct C = A - SIGMA*I in real arithmetic, and  |
c     | factor C in real arithmetic using LAPACK           |
c     | subroutine dgttrf. The matrix A is chosen to be    |
c     | the tridiagonal matrix derived from standard       |
c     | central difference of the 1-d convection diffusion |
c     | operator u" + rho*u' on the interval [0, 1] with   |
c     | zero Dirichlet boundary condition.                 |
c     %----------------------------------------------------%
c
      do i = 1,n
         c(i,i) = c(i,i) - sigmar
      end do
c
      call dgetrf(n, n, c, maxn, ipiv, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf in _NDRV2.'
         print*, ' '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in DNAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in DNAUPD to start the Arnoldi iteration. |
c     %-----------------------------------------------------%
c
      lworkl = 3*ncv**2+6*ncv
      tol    = zero
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed. Mode 3 of DNAUPD is used      |
c     | (IPARAM(7) = 3).  All these options can be        |
c     | changed by the user. For details see the          |
c     | documentation in DNAUPD.                          |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 3

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
c
 20   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dnaupd ( ido, bmat, n, which, nev, tol, resid,
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        info )
c
         if ( ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform  y <--- OP*x = inv[A-SIGMA*I]*x   |
c           | The user should supply his/her own linear |
c           | system solver here that takes             |
c           | workd(ipntr(1)) as the input, and returns |
c           | the result to workd(ipntr(2)).            |
c           %-------------------------------------------%
c
            call dcopy( n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
c
            call dgetrs('N', n, 1, c, maxn, ipiv,
     &                  workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV2.'
               print*, ' '
               go to 9000
            end if
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 20
c
         end if
c
c     %-----------------------------------------%
c     | Either we have convergence, or there is |
c     | an error.                               |
c     %-----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in DNAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation in _naupd.'
         print *, ' '
c
      else
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DNEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv,
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd,
     &        workl, lworkl, ierr )
c
c        %-----------------------------------------------%
c        | The real part of the eigenvalue is returned   |
c        | in the first column of the two dimensional    |
c        | array D, and the imaginary part is returned   |
c        | in the second column of D.  The corresponding |
c        | eigenvectors are returned in the first NEV    |
c        | columns of the two dimensional array V if     |
c        | requested.  Otherwise, an orthogonal basis    |
c        | for the invariant subspace corresponding to   |
c        | the eigenvalues in D is returned in V.        |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0 ) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of DNEUPD. |
c           %------------------------------------%
c
            print *, ' '
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
c
         else
c
            first  = .true.
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
               if (d(j,2) .eq. zero)  then
c
c                 %--------------------%
c                 | Ritz value is real |
c                 %--------------------%
c
                  call dgemv('N',n,n,1.0d0,a,maxn,v(1,j),1,0.0d0,ax,1)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
c
               else if (first) then
c
c                 %------------------------%
c                 | Ritz value is complex  |
c                 | Residual of one Ritz   |
c                 | value of the conjugate |
c                 | pair is computed.      |
c                 %------------------------%
c
                  call dgemv('N',n,n,1.0d0,a,maxn,v(1,j),1,0.0d0,ax,1)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                 call dgemv('N',n,n,1.0d0,a,maxn,v(1,j+1),1,0.0d0,ax,1)
                  call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
c
 30         continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call dmout(6, nconv, 3, d, maxncv, -6,
     &           'Ritz values (Real,Imag) and relative residuals')
c
         end if
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit',
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
c
         print *, ' '
         print *, ' _NDRV2 '
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
      if (isnan(v(1,1))) then
         stop 1
      end if
c
c     %---------------------------%
c     | Done with program dndrv2. |
c     %---------------------------%
c
 9000 continue
c
      end
