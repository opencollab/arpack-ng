      program cndrv1 
c
c     Example program to illustrate the idea of reverse communication
c     for a standard complex nonsymmetric eigenvalue problem. 
c
c     We implement example one of ex-complex.doc in DOCUMENTS directory
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is obtained from the standard central difference
c         discretization of the convection-diffusion operator 
c                 (Laplacian u) + rho*(du / dx)
c         on the unit squre [0,1]x[0,1] with zero Dirichlet boundary
c         condition.
c
c     ... OP = A  and  B = I.
c
c     ... Assume "call av (nx,x,y)" computes y = A*x
c
c     ... Use mode 1 of CNAUPD.
c
c\BeginLib
c
c\Routines called
c     cnaupd  ARPACK reverse communication interface routine.
c     cneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     scnrm2  Level 1 BLAS that computes the norm of a complex vector.
c     caxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     av      Matrix vector multiplication routine that computes A*x.
c     tv      Matrix vector multiplication routine that computes T*x,
c             where T is a tridiagonal matrix.  It is used in routine
c             av.
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
c FILE: ndrv1.F   SID: 2.4   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c---------------------------------------------------------------------------
c
c     %-----------------------------%
c     | Define maximum dimensions   |
c     | for all arrays.             |
c     | MAXN:   Maximum dimension   |
c     |         of the A allowed.   |
c     | MAXNEV: Maximum NEV allowed |
c     | MAXNCV: Maximum NCV allowed |
c     %-----------------------------%
c
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Complex 
     &                  ax(maxn), d(maxncv), 
     &                  v(ldv,maxncv), workd(3*maxn), 
     &                  workev(3*maxncv), resid(maxn), 
     &                  workl(3*maxncv*maxncv+5*maxncv)
      Real  
     &                  rwork(maxncv), rd(maxncv,3)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
      Complex 
     &                  sigma
      Real 
     &                  tol
      logical           rvec
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real 
     &                  scnrm2, slapy2
      external          scnrm2, caxpy, slapy2 
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %--------------------------------------------------%
c     | The number NX is the number of interior points   |
c     | in the discretization of the 2-dimensional       |
c     | convection-diffusion operator on the unit        |
c     | square with zero Dirichlet boundary condition.   | 
c     | The number N(=NX*NX) is the dimension of the     |
c     | matrix.  A standard eigenvalue problem is        |
c     | solved (BMAT = 'I').  NEV is the number of       |
c     | eigenvalues to be approximated.  The user can    |
c     | modify NX, NEV, NCV, WHICH to solve problems of  |
c     | different sizes, and to get different parts of   |
c     | the spectrum.  However, The following            |
c     | conditions must be satisfied:                    |
c     |                   N <= MAXN                      |
c     |                 NEV <= MAXNEV                    |
c     |           NEV + 2 <= NCV <= MAXNCV               | 
c     %--------------------------------------------------% 
c
      nx    = 10 
      n     = nx*nx 
      nev   = 4
      ncv   = 20 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'LM'
c
c     %---------------------------------------------------%
c     | The work array WORKL is used in CNAUPD as         | 
c     | workspace.  Its dimension LWORKL is set as        |
c     | illustrated below.  The parameter TOL determines  |
c     | the stopping criterion. If TOL<=0, machine        |
c     | precision is used.  The variable IDO is used for  |
c     | reverse communication, and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is  |
c     | generated to start the ARNOLDI iteration.         | 
c     %---------------------------------------------------%
c
      lworkl  = 3*ncv**2+5*ncv 
      tol    = 0.0 
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shift with respect to     |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of CNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | CNAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 1
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
c        | Repeatedly call the routine CNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call cnaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork,info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- OP*x                |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector, and return the matrix vector      |
c           | product to workd(ipntr(2)).               | 
c           %-------------------------------------------%
c
            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call CNAUPD again. |
c           %-----------------------------------------%
c
            go to 10

         end if
c 
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in CNAUPD  |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using CNEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call cneupd (rvec, 'A', select, d, v, ldv, sigma, 
     &        workev, bmat, n, which, nev, tol, resid, ncv, 
     &        v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &        rwork, ierr)
c
c        %----------------------------------------------%
c        | Eigenvalues are returned in the one          |
c        | dimensional array D.  The corresponding      |
c        | eigenvectors are returned in the first NCONV |
c        | (=IPARAM(5)) columns of the two dimensional  | 
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%
c
         if ( ierr .ne. 0) then
c 
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of CNEUPD. |
c           %------------------------------------%
c
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
c
         else
c
             nconv = iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
                call av(nx, v(1,j), ax)
                call caxpy(n, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = real (d(j))
                rd(j,2) = aimag(d(j))
                rd(j,3) = scnrm2(n, ax, 1)
                rd(j,3) = rd(j,3) / slapy2(rd(j,1),rd(j,2))
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call smout(6, nconv, 3, rd, maxncv, -6,
     &            'Ritz values (Real, Imag) and relative residuals')
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
         print *, '_NDRV1'
         print *, '====== '
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
c     | Done with program cndrv1. |
c     %---------------------------%
c
 9000 continue
c
      end
c 
c==========================================================================
c
c     matrix vector subroutine
c
c     The matrix used is the convection-diffusion operator
c     discretized using centered difference.
c
      subroutine av (nx, v, w)
      integer           nx, j, lo
      Complex          
     &                  v(nx*nx), w(nx*nx), one, h2
      parameter         (one = (1.0E+0, 0.0E+0) )
      external          caxpy, tv
c
c     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block 
c     tridiagonal matrix
c
c                  | T -I          | 
c                  |-I  T -I       |
c             OP = |   -I  T       |
c                  |        ...  -I|
c                  |           -I T|
c
c     derived from the standard central difference  discretization 
c     of the convection-diffusion operator (Laplacian u) + rho*(du/dx)
c     with zero boundary condition.
c
c     The subroutine TV is called to computed y<---T*x.
c
c
      h2 = one / cmplx((nx+1)*(nx+1))
c
      call tv(nx,v(1),w(1))
      call caxpy(nx, -one/h2, v(nx+1), 1, w(1), 1)
c
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call caxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
         call caxpy(nx, -one/h2, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
c
      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call caxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
c
      return
      end
c=========================================================================
      subroutine tv (nx, x, y)
c
      integer           nx, j 
      Complex 
     &                  x(nx), y(nx), h, h2, dd, dl, du
c
      Complex 
     &                  one, rho
      parameter         (one = (1.0E+0, 0.0E+0) ,
     &                   rho = (1.0E+2, 0.0E+0) )
c
c     Compute the matrix vector multiplication y<---T*x
c     where T is a nx by nx tridiagonal matrix with DD on the 
c     diagonal, DL on the subdiagonal, and DU on the superdiagonal
c     
      h   = one / cmplx(nx+1)
      h2  = h*h
      dd  = (4.0E+0, 0.0E+0)  / h2
      dl  = -one/h2 - (5.0E-1, 0.0E+0) *rho/h
      du  = -one/h2 + (5.0E-1, 0.0E+0) *rho/h
c 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end
