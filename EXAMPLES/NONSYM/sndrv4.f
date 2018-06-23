      program sndrv4
c
c     Simple program to illustrate the idea of reverse communication
c     in shift-invert mode for a generalized nonsymmetric eigenvalue
c     problem.
c
c     We implement example four of ex-nonsym.doc in DOCUMENTS directory
c
c\Example-4
c     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
c         where A and B are derived from the finite element discretization
c         of the 1-dimensional convection-diffusion operator
c                           (d^2u / dx^2) + rho*(du/dx)
c         on the interval [0,1] with zero Dirichlet boundary condition
c         using linear elements.
c
c     ... The shift sigma is a real number.
c
c     ... OP = inv[A-SIGMA*M]*M  and  B = M.
c
c     ... Use mode 3 of SNAUPD.
c
c\BeginLib
c
c\Routines called:
c     snaupd  ARPACK reverse communication interface routine.
c     sneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     sgttrf  LAPACK tridiagonal factorization routine.
c     sgttrs  LAPACK tridiagonal linear system solve routine.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     scopy   Level 1 BLAS that copies one vector to another.
c     sdot    Level 1 BLAS that computes the dot product of two vectors.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     av      Matrix vector multiplication routine that computes A*x.
c     mv      Matrix vector multiplication routine that computes M*x.
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
c FILE: ndrv4.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c-----------------------------------------------------------------------
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
      parameter         (maxn=256, maxnev=10, maxncv=25,
     &                   ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14), ipiv(maxn)
      logical           select(maxncv)
      Real
     &                  ax(maxn), mx(maxn), d(maxncv,3), resid(maxn),
     &                  v(ldv,maxncv), workd(3*maxn), workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv),
     &                  dd(maxn), dl(maxn), du(maxn),
     &                  du2(maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, ierr, j,
     &                  nconv, maxitr, ishfts, mode
      Real
     &                  tol, h, s,
     &                  sigmar, sigmai, s1, s2, s3
      logical           first, rvec
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real
     &                  sdot, snrm2, slapy2
      external          sdot, snrm2, slapy2, sgttrf, sgttrs
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic         abs
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &                   one, zero, two, six, rho
      common            /convct/ rho
      parameter         (one = 1.0E+0, zero = 0.0E+0,
     &                   two = 2.0E+0, six = 6.0E+0)
c
c     %-----------------------%
c     | Executable statements |
c     %-----------------------%
c
c     %----------------------------------------------------%
c     | The number N is the dimension of the matrix.  A    |
c     | generalized eigenvalue problem is solved (BMAT =   |
c     | 'G').  NEV is the number of eigenvalues (closest   |
c     | to SIGMAR) to be approximated.  Since the          |
c     | shift-invert mode is used,  WHICH is set to 'LM'.  |
c     | The user can modify NEV, NCV, SIGMAR to solve      |
c     | problems of different sizes, and to get different  |
c     | parts of the spectrum.  However, The following     |
c     | conditions must be satisfied:                      |
c     |                     N <= MAXN,                     |
c     |                   NEV <= MAXNEV,                   |
c     |               NEV + 2 <= NCV <= MAXNCV             |
c     %----------------------------------------------------%
c
      n     = 100
      nev   = 4
      ncv   = 10
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV4: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV4: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV4: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
      sigmar = one
      sigmai = zero
c
c     %--------------------------------------------------%
c     | Construct C = A - SIGMA*M in real arithmetic,    |
c     | and factor C in real arithmetic (using LAPACK    |
c     | subroutine sgttrf). The matrix A is chosen to be |
c     | the tridiagonal matrix derived from the standard |
c     | central difference discretization of the 1-d     |
c     | convection-diffusion operator u" + rho*u' on the |
c     | interval [0, 1] with zero Dirichlet boundary     |
c     | condition.  The matrix M is the mass matrix      |
c     | formed by using piecewise linear elements on     |
c     | [0,1].                                           |
c     %--------------------------------------------------%
c
      rho = 1.0E+1
      h = one / real(n+1)
      s = rho / two
c
      s1 = -one/h - s - sigmar*h/six
      s2 = two/h - 4.0E+0*sigmar*h/six
      s3 = -one/h + s - sigmar*h/six
c
      do 10 j = 1, n-1
         dl(j) = s1
         dd(j) = s2
         du(j) = s3
  10  continue
      dd(n) = s2
c
      call sgttrf(n, dl, dd, du, du2, ipiv, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf in _NDRV4.'
         print*, ' '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in SNAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in SNAUPD to start the Arnoldi iteration. |
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
c     | iterations allowed.  Mode 3 of SNAUPD is used     |
c     | (IPARAM(7) = 3).  All these options can be        |
c     | changed by the user. For details, see the         |
c     | documentation in SNAUPD.                          |
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
c     %------------------------------------------%
c     | M A I N   L O O P(Reverse communication) |
c     %------------------------------------------%
c
 20   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine SNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call snaupd ( ido, bmat, n, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd,
     &                 workl, lworkl, info )
c
         if (ido .eq. -1) then
c
c           %-------------------------------------------%
c           | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x |
c           | to force starting vector into the range   |
c           | of OP.   The user should supply his/her   |
c           | own matrix vector multiplication routine  |
c           | and a linear system solver.  The matrix   |
c           | vector multiplication routine should take |
c           | workd(ipntr(1)) as the input. The final   |
c           | result should be returned to              |
c           | workd(ipntr(2)).                          |
c           %-------------------------------------------%
c
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
            call sgttrs('N', n, 1, dl, dd, du, du2, ipiv,
     &                  workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV4.'
               print*, ' '
               go to 9000
            end if
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SNAUPD again. |
c           %-----------------------------------------%
c
            go to 20
c
         else if ( ido .eq. 1) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x = inv[A-sigma*M]*M*x |
c           | M*x has been saved in workd(ipntr(3)).  |
c           | The user only need the linear system    |
c           | solver here that takes workd(ipntr(3))  |
c           | as input, and returns the result to     |
c           | workd(ipntr(2)).                        |
c           %-----------------------------------------%
c
            call scopy( n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
            call sgttrs ('N', n, 1, dl, dd, du, du2, ipiv,
     &                  workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV4.'
               print*, ' '
               go to 9000
            end if
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SNAUPD again. |
c           %-----------------------------------------%
c
            go to 20
c
         else if ( ido .eq. 2) then
c
c           %---------------------------------------------%
c           |          Perform  y <--- M*x                |
c           | Need matrix vector multiplication routine   |
c           | here that takes workd(ipntr(1)) as input    |
c           | and returns the result to workd(ipntr(2)).  |
c           %---------------------------------------------%
c
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SNAUPD again. |
c           %-----------------------------------------%
c
            go to 20
c
         end if
c
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
c        | documentation in SNAUPD. |
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
c        | Post-Process using SNEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
         call sneupd ( rvec, 'A', select, d, d(1,2), v, ldv,
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd,
     &        workl, lworkl, ierr )
c
c        %-----------------------------------------------%
c        | The real part of the eigenvalue is returned   |
c        | in the first column of the two dimensional    |
c        | array D, and the IMAGINARY part is returned   |
c        | in the second column of D.  The corresponding |
c        | eigenvectors are returned in the first NEV    |
c        | columns of the two dimensional array V if     |
c        | requested.  Otherwise, an orthogonal basis    |
c        | for the invariant subspace corresponding to   |
c        | the eigenvalues in D is returned in V.        |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of SNEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
c
         else
c
             first = .true.
             nconv =  iparam(5)
             do 30 j=1, nconv
c
c                %---------------------------%
c                | Compute the residual norm |
c                |                           |
c                |   ||  A*x - lambda*x ||   |
c                |                           |
c                | for the NCONV accurately  |
c                | computed eigenvalues and  |
c                | eigenvectors.  (iparam(5) |
c                | indicates how many are    |
c                | accurate to the requested |
c                | tolerance)                |
c                %---------------------------%
c
                 if (d(j,2) .eq. zero)  then
c
c                    %--------------------%
c                    | Ritz value is real |
c                    %--------------------%
c
                     call av(n, v(1,j), ax)
                     call mv(n, v(1,j), mx)
                     call saxpy(n, -d(j,1), mx, 1, ax, 1)
                     d(j,3) = snrm2(n, ax, 1)
                     d(j,3) = d(j,3) / abs(d(j,1))
c
                 else if (first) then
c
c                    %------------------------%
c                    | Ritz value is complex. |
c                    | Residual of one Ritz   |
c                    | value of the conjugate |
c                    | pair is computed.      |
c                    %------------------------%
c
                     call av(n, v(1,j), ax)
                     call mv(n, v(1,j), mx)
                     call saxpy(n, -d(j,1), mx, 1, ax, 1)
                     call mv(n, v(1,j+1), mx)
                     call saxpy(n, d(j,2), mx, 1, ax, 1)
                     d(j,3) = snrm2(n, ax, 1)
                     call av(n, v(1,j+1), ax)
                     call mv(n, v(1,j+1), mx)
                     call saxpy(n, -d(j,1), mx, 1, ax, 1)
                     call mv(n, v(1,j), mx)
                     call saxpy(n, -d(j,2), mx, 1, ax, 1)
                     d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
                     d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
                     d(j+1,3) = d(j,3)
                     first = .false.
                 else
                     first = .true.
                 end if
c
  30         continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call smout(6, nconv, 3, d, maxncv, -6,
     &            'Ritz values (Real,Imag) and relative residuals')
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
         print *, ' _NDRV4 '
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
c     | Done with program sndrv4. |
c     %---------------------------%
c
 9000 continue
c
      end
c
c==========================================================================
c
c     matrix vector multiplication subroutine
c
      subroutine mv (n, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), one, four, six, h
      parameter         (one = 1.0E+0, four = 4.0E+0, six = 6.0E+0)
c
c     Compute the matrix vector multiplication y<---M*x
c     where M is mass matrix formed by using piecewise linear elements
c     on [0,1].
c
      w(1) =  ( four*v(1) + one*v(2) ) / six
      do 10 j = 2,n-1
         w(j) = ( one*v(j-1) + four*v(j) + one*v(j+1) ) / six
 10   continue
      w(n) =  ( one*v(n-1) + four*v(n) ) / six
c
      h = one / real(n+1)
      call sscal(n, h, w, 1)
      return
      end
c------------------------------------------------------------------
      subroutine av (n, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), one, two, dd, dl, du, s, h, rho
      common            /convct/ rho
      parameter         (one = 1.0E+0, two = 2.0E+0)
c
c     Compute the matrix vector multiplication y<---A*x
c     where A is obtained from the finite element discretization of the
c     1-dimensional convection diffusion operator
c                     d^u/dx^2 + rho*(du/dx)
c     on the interval [0,1] with zero Dirichlet boundary condition
c     using linear elements.
c     This routine is only used in residual calculation.
c
      h = one / real(n+1)
      s = rho / two
      dd = two / h
      dl = -one/h - s
      du = -one/h + s
c
      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(n) =  dl*v(n-1) + dd*v(n)
      return
      end

