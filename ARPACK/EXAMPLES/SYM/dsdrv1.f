      program dsdrv1 
c
c     Simple program to illustrate the idea of reverse communication
c     in regular mode for a standard symmetric eigenvalue problem.
c
c     We implement example one of ex-sym.doc in SRC directory
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is derived from the central difference discretization
c         of the 2-dimensional Laplacian on the unit square [0,1]x[0,1]
c         with zero Dirichlet boundary condition.
c
c     ... OP = A  and  B = I.
c
c     ... Assume "call av (n,x,y)" computes y = A*x.
c
c     ... Use mode 1 of DSAUPD.
c
c\BeginLib
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
c FILE: sdrv1.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
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
      integer          maxn, maxnev, maxncv, ldv
      parameter       (maxn=256, maxnev=10, maxncv=25, 
     $                 ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Double precision
     &                 v(ldv,maxncv), workl(maxncv*(maxncv+8)),
     &                 workd(3*maxn), d(maxncv,2), resid(maxn),
     &                 ax(maxn)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr, j, 
     &                 nx, nconv, maxitr, mode, ishfts
      logical          rvec
      Double precision      
     &                 tol, sigma
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                 zero
      parameter        (zero = 0.0D+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision           
     &                 dnrm2
      external         dnrm2, daxpy
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
c     | The number NX is the number of interior points     |
c     | in the discretization of the 2-dimensional         |
c     | Laplacian on the unit square with zero Dirichlet   |
c     | boundary condition.  The number N(=NX*NX) is the   |
c     | dimension of the matrix.  A standard eigenvalue    |
c     | problem is solved (BMAT = 'I'). NEV is the number  |
c     | of eigenvalues to be approximated.  The user can   |
c     | modify NEV, NCV, WHICH to solve problems of        |
c     | different sizes, and to get different parts of the |
c     | spectrum.  However, The following conditions must  |
c     | be satisfied:                                      |
c     |                   N <= MAXN,                       | 
c     |                 NEV <= MAXNEV,                     |
c     |             NEV + 1 <= NCV <= MAXNCV               | 
c     %----------------------------------------------------% 
c
      nx = 10
      n = nx*nx
      nev =  4 
      ncv =  10 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'SM'
c
c     %--------------------------------------------------%
c     | The work array WORKL is used in DSAUPD as        |
c     | workspace.  Its dimension LWORKL is set as       |
c     | illustrated below.  The parameter TOL determines |
c     | the stopping criterion.  If TOL<=0, machine      |
c     | precision is used.  The variable IDO is used for |
c     | reverse communication and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is |
c     | generated in DSAUPD to start the Arnoldi         |
c     | iteration.                                       |
c     %--------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol = zero 
      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1).  All these options may be        |
c     | changed by the user. For details, see the         |
c     | documentation in DSAUPD.                          |
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
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, workl,
     &                 lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c
            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
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
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        %-------------------------------------------%
c           
         rvec = .true.
c
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )

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
c
         if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of DSEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
c
         else
c
             nconv =  iparam(5)
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
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
c
 20          continue
c
c            %-------------------------------%
c            | Display computed residuals    |
c            %-------------------------------%
c
             call dmout(6, nconv, 2, d, maxncv, -6,
     &            'Ritz values and relative residuals')
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
         print *, ' _SDRV1 '
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
c     | Done with program dsdrv1. |
c     %---------------------------%
c
 9000 continue
c
      end
c 
c ------------------------------------------------------------------
c     matrix vector subroutine
c
c     The matrix used is the 2 dimensional discrete Laplacian on unit
c     square with zero Dirichlet boundary condition.
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
c     The subroutine TV is called to computed y<---T*x.
c
      subroutine av (nx, v, w)
      integer           nx, j, lo, n2
      Double precision
     &                  v(nx*nx), w(nx*nx), one, h2
      parameter         ( one = 1.0D+0 ) 
c
      call tv(nx,v(1),w(1))
      call daxpy(nx, -one, v(nx+1), 1, w(1), 1)
c
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
         call daxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
c
      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call daxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
c
c     Scale the vector w by (1/h^2), where h is the mesh size
c
      n2 = nx*nx
      h2 = one / dble((nx+1)*(nx+1))
      call dscal(n2, one/h2, w, 1) 
      return
      end
c
c-------------------------------------------------------------------
      subroutine tv (nx, x, y)
c
      integer           nx, j 
      Double precision
     &                  x(nx), y(nx), dd, dl, du
c
      Double precision
     &                 one
      parameter        (one = 1.0D+0 )
c
c     Compute the matrix vector multiplication y<---T*x
c     where T is a nx by nx tridiagonal matrix with DD on the 
c     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
c     
c
      dd  = 4.0D+0
      dl  = -one 
      du  = -one
c 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end

