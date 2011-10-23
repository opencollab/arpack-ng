      program sndrv3
c
c     Simple program to illustrate the idea of reverse communication
c     in inverse mode for a generalized nonsymmetric eigenvalue problem.
c
c     We implement example three of ex-nonsym.doc in DOCUMENTS directory
c
c\Example-3
c     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
c         where A and B are derived from the finite element discretization
c         of the 1-dimensional convection-diffusion operator
c                           (d^2u / dx^2) + rho*(du/dx)
c         on the interval [0,1] with zero Dirichlet boundary condition
c         using linear elements.
c
c     ... So OP = inv[M]*A  and  B = M.
c
c     ... Use mode 2 of SNAUPD. 
c
c\BeginLib
c
c\Routines called:
c     snaupd  ARPACK reverse communication interface routine.
c     sneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     spttrf  LAPACK symmetric positive definite tridiagonal factorization
c             routine.
c     spttrs  LAPACK symmetric positive definite tridiagonal solve routine.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
c FILE: ndrv3.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
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
      parameter         (maxn=256, maxnev=10, maxncv=25, 
     &                   ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Real
     &                  ax(maxn), mx(maxn), d(maxncv, 3), resid(maxn),
     &                  v(ldv,maxncv), workd(3*maxn), 
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv),
     &                  md(maxn), me(maxn-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, ierr, j,
     &                  nconv, maxitr, ishfts, mode
      Real           
     &                  tol, sigmar, sigmai, h
      logical           first, rvec
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real  
     &                  zero, one
      parameter         (zero = 0.0E+0, one = 1.0E+0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
      Real 
     &                  snrm2, slapy2
      external          saxpy, snrm2, spttrf, spttrs, slapy2
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic         abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %----------------------------------------------------%
c     | The number N is the dimension of the matrix.  A    |
c     | generalized eigenvalue problem is solved (BMAT =   |
c     | 'G').  NEV is the number of eigenvalues to be      |
c     | approximated.  The user can modify NEV, NCV, WHICH |
c     | to solve problems of different sizes, and to get   |
c     | different parts of the spectrum.  However, The     |
c     | following conditions must be satisfied:            |
c     |                    N <= MAXN,                      |
c     |                  NEV <= MAXNEV,                    |
c     |              NEV + 2 <= NCV <= MAXNCV              |
c     %----------------------------------------------------%
c
      n     = 100 
      nev   = 4 
      ncv   = 20 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV3: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
c
c     %------------------------------------------------%
c     | M is the mass matrix formed by using piecewise |
c     | linear elements on [0,1].                      |
c     %------------------------------------------------%
c
      h = one / real(n+1)
      do 20 j = 1, n-1
         md(j) = 4.0E+0*h 
         me(j) = one*h 
  20  continue 
      md(n) = 4.0E+0*h
c 
      call spttrf(n, md, me, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _pttrf. ' 
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
      tol    = 0.0 
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 2 of SNAUPD is used     |
c     | (IPARAM(7) = 2).  All these options can be        |
c     | changed by the user. For details, see the         |
c     | documentation in SNAUPD.                          |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 2
c
      iparam(1) = ishfts 
      iparam(3) = maxitr  
      iparam(7) = mode 
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
c
 10   continue
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
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %----------------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x      |
c           | The user should supply his/her own     |
c           | matrix vector routine and a linear     |
c           | system solver.  The matrix-vector      |
c           | subroutine should take workd(ipntr(1)) |
c           | as input, and the final result should  |
c           | be returned to workd(ipntr(2)).        |
c           %----------------------------------------%
c
            call av (n, workd(ipntr(1)), workd(ipntr(2)))
            call spttrs(n, 1, md, me, workd(ipntr(2)), n,
     &                  ierr)
            if ( ierr .ne. 0 ) then
               print*, ' ' 
               print*, ' ERROR with _pttrs. ' 
               print*, ' '
               go to 9000
            end if
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         else if ( ido .eq. 2) then
c
c           %-------------------------------------%
c           |        Perform  y <--- M*x          |
c           | The matrix vector multiplication    |
c           | routine should take workd(ipntr(1)) |
c           | as input and return the result to   |
c           | workd(ipntr(2)).                    |
c           %-------------------------------------%
c
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call SNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
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
c        | Error message. Check the |
c        | documentation in SNAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd.'
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
         if ( ierr .ne. 0 ) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of SNEUPD. |
c           %------------------------------------%
c 
            print *, ' ' 
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd'
            print *, ' '
c
         else 
c
            first = .true.
            nconv = iparam(5)
            do 30 j=1, iparam(5)
c
c              %---------------------------%
c              | Compute the residual norm |
c              |                           |
c              |  ||  A*x - lambda*M*x ||  |
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
                  call av(n, v(1,j), ax)
                  call mv(n, v(1,j), mx)
                  call saxpy(n, -d(j,1), mx, 1, ax, 1)
                  d(j,3) = snrm2(n, ax, 1)
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
                  call av(n, v(1,j), ax)
                  call mv(n, v(1,j), mx)
                  call saxpy(n, -d(j,1), mx, 1, ax, 1)
                  call mv(n, v(1,j+1), mx)
                  call saxpy(n, d(j,2), mx, 1, ax, 1)
                  d(j,3) = snrm2(n, ax, 1)**2
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
  30        continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call smout(6, nconv, 3, d, maxncv, -6,
     &           'Ritz values (Real,Imag) and relative residuals')
c
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
         print *, ' _NDRV3 '
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
c     | Done with program sndrv3. |
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
      subroutine av (n, v, w)
      integer           n, j
      Real            
     &                  v(n), w(n), one, two, dd, dl, du, s, h, rho 
      parameter         ( rho = 1.0E+1, one = 1.0E+0, 
     &                    two = 2.0E+0)
c
c     Compute the matrix vector multiplication y<---A*x
c     where A is stiffness matrix obtained from the finite element
c     discretization of the 1-dimensional convection diffusion operator
c                           d^2u/dx^2 + rho*(du/dx)
c     on the interval [0,1] with zero Dirichlet boundary condition using
c     linear elements.
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
c------------------------------------------------------------------------
      subroutine mv (n, v, w)
      integer           n, j
      Real 
     &                  v(n), w(n), one, four, h
      parameter         ( one = 1.0E+0, four = 4.0E+0)
c
c     Compute the matrix vector multiplication y<---M*x
c     where M is the mass matrix formed by using piecewise linear 
c     elements on [0,1].
c 
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,n-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1) 
 10   continue 
      w(n) =  one*v(n-1) + four*v(n) 
c
      h = one / real(n+1)
      call sscal(n, h, w, 1)
      return
      end
