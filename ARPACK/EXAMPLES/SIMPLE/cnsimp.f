      program cnsimp 
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
c     where A is a general n by n complex matrix.
c
c     The main points illustrated here are 
c
c        1) How to declare sufficient memory to find NEV 
c           eigenvalues of largest magnitude.  Other options
c           are available.
c
c        2) Illustration of the reverse communication interface 
c           needed to utilize the top level ARPACK routine CNAUPD 
c           that computes the quantities needed to construct
c           the desired eigenvalues and eigenvectors(if requested).
c
c        3) How to extract the desired eigenvalues and eigenvectors
c           using the ARPACK routine CNEUPD.
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
c
c     Once usage of this routine is understood, you may wish to explore
c     the other available options to improve convergence, to solve generalized
c     problems, etc.  Look at the file ex-complex.doc in DOCUMENTS directory.
c     This codes implements  
c
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c     ... OP = A  and  B = I.
c     ... Assume "call av (nx,x,y)" computes y = A*x
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
c FILE: nsimp.F   SID: 2.4   DATE OF SID: 10/20/00   RELEASE: 2
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
      Complex 
     &                  ax(maxn), d(maxncv), 
     &                  v(ldv,maxncv), workd(3*maxn), 
     &                  workev(2*maxncv), resid(maxn), 
     &                  workl(3*maxncv*maxncv+5*maxncv)
      Real  
     &                  rwork(maxncv), rd(maxncv,3)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, ierr,
     &                  j, ishfts, maxitr, mode1, nconv
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
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting mcaupd = 1                     |
c     %-------------------------------------------------%
c
      include 'debug.h'
      ndigit = -3
      logfil = 6
      mcaitr = 0 
      mcapps = 0
      mcaupd = 1
      mcaup2 = 0
      mceigh = 0
      mceupd = 0
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
c     |       factorization                           |
c     |                                               |
c     |    3) This is a standard problem              |
c     |         (indicated by bmat  = 'I')            |
c     |                                               |
c     |    4) Ask for the NEV eigenvalues of          |
c     |       largest magnitude                       |
c     |         (indicated by which = 'LM')           |
c     |       See documentation in CNAUPD for the     |
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
c     | conditions before calling CNAUPD                    |
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
c     |      from CNAUPD. (see usage below)                 |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to CNAUPD.                                | 
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
c     | The work array WORKL is used in CNAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl  = 3*ncv**2+5*ncv 
      tol    = 0.0 
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting IPARAM(1) = 1).             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of CNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | CNAUPD.                                           |
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
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse Communication Loop) | 
c     %------------------------------------------------%
c
 10   continue
c   
c        %---------------------------------------------%
c        | Repeatedly call the routine CNAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%

         call cnaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &                 v, ldv, iparam, ipntr, workd, workl, lworkl,
     &                 rwork,info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                                           |
c           |                y <--- A*x                 |
c           |                                           |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector x , and returns the resulting      |
c           | matrix-vector product y = A*x in the      |
c           | array workd(ipntr(2)).                    | 
c           %-------------------------------------------%
c
            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
c 
c           %-----------------------------------------%
c           | L O O P   B A C K to call CNAUPD again. |
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
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine CNEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c           
         rvec = .true.
c
         call cneupd (rvec, 'A', select, D, V, ldv, sigma,
     &        workev, bmat, n, which, nev, tol, resid, ncv,
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,
     &        rwork, ierr)
c
c        %-----------------------------------------------%
c        | Eigenvalues are returned in the one           |
c        | dimensional array D and the corresponding     |
c        | eigenvectors are returned in the first        |
c        | NCONV (=IPARAM(5)) columns of the two         |
c        | dimensional array V if requested.  Otherwise, |
c        | an orthogonal basis for the invariant         |
c        | subspace corresponding to the eigenvalues in  |
c        | D is returned in V.                           |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of CNEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
c
         else
c
             nconv =  iparam(5)
             do 20 j=1, nconv
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
         print *, '_NSIMP '
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
c     | Done with program cnsimp. |
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
      external          caxpy
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
c     of the 2-dimensional convection-diffusion operator 
c                  (Laplacian u) + rho*(du/dx)
c     on the unit squqre with zero boundary condition.
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
