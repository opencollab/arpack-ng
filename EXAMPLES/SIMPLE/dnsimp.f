      program dnsimp
c
c
c     This example program is intended to illustrate the
c     simplest case of using ARPACK in considerable detail.
c     This code may be used to understand basic usage of ARPACK
c     and as a template for creating an interface to ARPACK.
c
c     This code shows how to use ARPACK to find a few eigenvalues
c     (lambda) and corresponding eigenvectors (x) for the standard
c     eigenvalue problem:
c
c                        A*x = lambda*x
c
c     where A is a n by n real nonsymmetric matrix.
c
c     The main points illustrated here are
c
c        1) How to declare sufficient memory to find NEV
c           eigenvalues of largest magnitude.  Other options
c           are available.
c
c        2) Illustration of the reverse communication interface
c           needed to utilize the top level ARPACK routine DNAUPD
c           that computes the quantities needed to construct
c           the desired eigenvalues and eigenvectors(if requested).
c
c        3) How to extract the desired eigenvalues and eigenvectors
c           using the ARPACK routine DNEUPD.
c
c     The only thing that must be supplied in order to use this
c     routine on your problem is to change the array dimensions
c     appropriately, to specify WHICH eigenvalues you want to compute
c     and to supply a matrix-vector product
c
c                         w <-  Av
c
c     in place of the call to AV( )  below.
c
c     Once usage of this routine is understood, you may wish to explore
c     the other available options to improve convergence, to solve generalized
c     problems, etc.  Look at the file ex-nonsym.doc in DOCUMENTS directory.
c     This codes implements
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is obtained from the standard central difference
c         discretization of the convection-diffusion operator 
c                 (Laplacian u) + rho*(du / dx)
c         on the unit square, with zero Dirichlet boundary condition.
c
c     ... OP = A  and  B = I.
c     ... Assume "call av (nx,x,y)" computes y = A*x
c     ... Use mode 1 of DNAUPD.
c
c\BeginLib
c
c\Routines called:
c     dnaupd  ARPACK reverse communication interface routine.
c     dneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
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
c FILE: nsimp.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c---------------------------------------------------------------------------
c
c     %------------------------------------------------------%
c     | Storage Declarations:                                |
c     |                                                      |
c     | The maximum dimensions for all arrays are            |
c     | set here to accommodate a problem size of            |
c     | N .le. MAXN                                          |
c     |                                                      |
c     | NEV is the number of eigenvalues requested.          |
c     |     See specifications for ARPACK usage below.       |
c     |                                                      |
c     | NCV is the largest number of basis vectors that will |
c     |     be used in the Implicitly Restarted Arnoldi      |
c     |     Process.  Work per major iteration is            |
c     |     proportional to N*NCV*NCV.                       |
c     |                                                      |
c     | You must set:                                        |
c     |                                                      |
c     | MAXN:   Maximum dimension of the A allowed.          |
c     | MAXNEV: Maximum NEV allowed.                         |
c     | MAXNCV: Maximum NCV allowed.                         |
c     %------------------------------------------------------%
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
      Double precision
     &                  ax(maxn), d(maxncv,3), resid(maxn), 
     &                  v(ldv,maxncv), workd(3*maxn), 
     &                  workev(3*maxncv), 
     &                  workl(3*maxncv*maxncv+6*maxncv)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, ierr,
     &                  j, ishfts, maxitr, mode1, nconv
      Double precision
     &                  tol, sigmar, sigmai
      logical           first, rvec
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                  zero
      parameter         (zero = 0.0D+0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision
     &                  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy 
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
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting mnaupd = 1.                    |
c     %-------------------------------------------------%
c
      include 'debug.h'
      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 1
      mnaup2 = 0
      mneigh = 0
      mneupd = 0
c
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
      nx    = 10 
      n     = nx*nx 
c
c     %-----------------------------------------------%
c     |                                               |
c     | Specifications for ARPACK usage are set       |
c     | below:                                        |
c     |                                               |
c     |    1) NEV = 4  asks for 4 eigenvalues to be   |
c     |       computed.                               |
c     |                                               |
c     |    2) NCV = 20 sets the length of the Arnoldi |
c     |       factorization.                          |
c     |                                               |
c     |    3) This is a standard problem.             |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude.                      |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in DNAUPD for the     |
c     |       other options SM, LR, SR, LI, SI.       |
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 2 <= NCV <= MAXNCV             |
c     |                                               |
c     %-----------------------------------------------%
c
      nev   = 4
      ncv   = 20
      bmat  = 'I'
      which = 'LM'
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling DNAUPD                    |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DNAUPD. (see usage below)                 |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DNAUPD.                                |
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     |
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID).  |
c     |                                                     |
c     | The work array WORKL is used in DNAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl  = 3*ncv**2+6*ncv 
      tol    = zero 
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1).             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DNAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode1 = 1
c
      iparam(1) = ishfts
c
      iparam(3) = maxitr
c
      iparam(7) = mode1
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) | 
c     %-------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, 
     &                 v, ldv, iparam, ipntr, workd, workl, lworkl, 
     &                 info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- Op*x                |
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
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         endif
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
c        | documentation in DNAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _naupd, info = ',info
         print *, ' Check the documentation of _naupd'
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
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        |                                           |
c        | The routine DNEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1,)                                   |
c        |                                           |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv, 
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
c
c        %------------------------------------------------%
c        | The real parts of the eigenvalues are returned |
c        | in the first column of the two dimensional     |
c        | array D, and the IMAGINARY part are returned   |
c        | in the second column of D.  The corresponding  |
c        | eigenvectors are returned in the first         |
c        | NCONV (= IPARAM(5)) columns of the two         |
c        | dimensional array V if requested.  Otherwise,  |
c        | an orthogonal basis for the invariant subspace |
c        | corresponding to the eigenvalues in D is       |
c        | returned in V.                                 |
c        %------------------------------------------------%
c
         if ( ierr .ne. 0) then
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
            first = .true.
            nconv =  iparam(5)
            do 20 j=1, nconv
c
c              %---------------------------%
c              | Compute the residual norm |
c              |                           |
c              |   ||  A*x - lambda*x ||   |
c              |                           |
c              | for the NCONV accurately  |
c              | computed eigenvalues and  |
c              | eigenvectors.  (IPARAM(5) |
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
                  call av(nx, v(1,j), ax)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
c
               else if (first) then
c
c                 %------------------------%
c                 | Ritz value is complex. |
c                 | Residual of one Ritz   |
c                 | value of the conjugate |
c                 | pair is computed.      |
c                 %------------------------%
c
                  call av(nx, v(1,j), ax)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  call av(nx, v(1,j+1), ax)
                  call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
c
 20         continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call dmout(6, nconv, 3, d, maxncv, -6,
     &           'Ritz values (Real, Imag) and residual residuals')
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
         print *, ' _NSIMP '
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
c     | Done with program dnsimp. |
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
c     The matrix used is the 2 dimensional convection-diffusion 
c     operator discretized using central difference.
c
      subroutine av (nx, v, w)
      integer           nx, j, lo
      Double precision         
     &                  v(nx*nx), w(nx*nx), one, h2
      parameter         (one = 1.0D+0)
      external          daxpy, tv
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
c     derived from the standard central difference discretization 
c     of the 2 dimensional convection-diffusion operator 
c     (Laplacian u) + rho*(du/dx) on a unit square with zero boundary 
c     condition.
c
c     When rho*h/2 <= 1, the discrete convection-diffusion operator 
c     has real eigenvalues.  When rho*h/2 > 1, it has complex 
c     eigenvalues.
c
c     The subroutine TV is called to computed y<---T*x.
c
c
      h2 = one / dble((nx+1)*(nx+1))
c
      call tv(nx,v(1),w(1))
      call daxpy(nx, -one/h2, v(nx+1), 1, w(1), 1)
c
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call daxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
         call daxpy(nx, -one/h2, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
c
      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call daxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
c
      return
      end
c=========================================================================
      subroutine tv (nx, x, y)
c
      integer           nx, j 
      Double precision
     &                  x(nx), y(nx), h, dd, dl, du, h2
c
      Double precision
     &                  one, rho
      parameter         (one = 1.0D+0, rho = 1.0D+2)
c
c     Compute the matrix vector multiplication y<---T*x
c     where T is a nx by nx tridiagonal matrix with DD on the 
c     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
c
c     When rho*h/2 <= 1, the discrete convection-diffusion operator 
c     has real eigenvalues.  When rho*h/2 > 1, it has complex 
c     eigenvalues.
c
      h   = one / dble(nx+1)
      h2  = h*h
      dd  = 4.0D+0 / h2 
      dl  = -one/h2 - 5.0D-1*rho/h
      du  = -one/h2 + 5.0D-1*rho/h
c 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end
