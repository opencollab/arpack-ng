      program sssimp
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
c     where A is an n by n real symmetric matrix.
c
c     The main points illustrated here are 
c
c        1) How to declare sufficient memory to find NEV 
c           eigenvalues of largest magnitude.  Other options
c           are available.
c
c        2) Illustration of the reverse communication interface 
c           needed to utilize the top level ARPACK routine SSAUPD 
c           that computes the quantities needed to construct
c           the desired eigenvalues and eigenvectors(if requested).
c
c        3) How to extract the desired eigenvalues and eigenvectors
c           using the ARPACK routine SSEUPD.
c
c     The only thing that must be supplied in order to use this
c     routine on your problem is to change the array dimensions 
c     appropriately, to specify WHICH eigenvalues you want to compute 
c     and to supply a matrix-vector product
c
c                         w <-  Av
c
c     in place of the call to AV( ) below.
c
c     Once usage of this routine is understood, you may wish to explore
c     the other available options to improve convergence, to solve generalized
c     problems, etc.  Look at the file ex-sym.doc in DOCUMENTS directory.
c     This codes implements  
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is derived from the central difference discretization
c         of the 2-dimensional Laplacian on the unit square with
c         zero Dirichlet boundary condition.
c     ... OP = A  and  B = I.
c     ... Assume "call av (n,x,y)" computes y = A*x
c     ... Use mode 1 of SSAUPD.
c
c\BeginLib
c
c\Routines called:
c     ssaupd  ARPACK reverse communication interface routine.
c     sseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
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
c FILE: ssimp.F   SID: 2.6   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
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
      integer          maxn, maxnev, maxncv, ldv
      parameter       (maxn=256, maxnev=10, maxncv=25, 
     $                 ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Real
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
      integer          ido, n, nev, ncv, lworkl, info, ierr,
     &                 j, nx, ishfts, maxitr, mode1, nconv
      logical          rvec
      Real      
     &                 tol, sigma
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &                 zero
      parameter        (zero = 0.0E+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real           
     &                 snrm2
      external         snrm2, saxpy
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
c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting msaupd = 1.                    |
c     %-------------------------------------------------%
c
      include 'debug.h'
      ndigit = -3
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      msaupd = 1
      msaup2 = 0
      mseigt = 0
      mseupd = 0
c     
c     %-------------------------------------------------%
c     | The following sets dimensions for this problem. |
c     %-------------------------------------------------%
c
      nx = 10
      n = nx*nx
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
c     |       See documentation in SSAUPD for the     |
c     |       other options SM, LA, SA, LI, SI.       | 
c     |                                               |
c     | Note: NEV and NCV must satisfy the following  |
c     | conditions:                                   |
c     |              NEV <= MAXNEV                    |
c     |          NEV + 1 <= NCV <= MAXNCV             |
c     %-----------------------------------------------%
c
      nev   = 4
      ncv   = 20 
      bmat  = 'I'
      which = 'LM'
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling SSAUPD                    |
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
c     |      from SSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to SSAUPD.                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  | 
c     |                                                     |
c     | The work array WORKL is used in SSAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol = zero 
      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting PARAM(1) = 1).              |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of SSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | SSAUPD.                                           |
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
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
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
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in SSAUPD. |
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
c        | Post-Process using SSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine SSEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c           
          rvec = .true.
c
          call sseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     &         bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &         iparam, ipntr, workd, workl, lworkl, ierr )
c
c         %----------------------------------------------%
c         | Eigenvalues are returned in the first column |
c         | of the two dimensional array D and the       |
c         | corresponding eigenvectors are returned in   |
c         | the first NCONV (=IPARAM(5)) columns of the  |
c         | two dimensional array V if requested.        |
c         | Otherwise, an orthogonal basis for the       |
c         | invariant subspace corresponding to the      |
c         | eigenvalues in D is returned in V.           |
c         %----------------------------------------------%
c
          if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of SSEUPD. |
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
                call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = snrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
c
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call smout(6, nconv, 2, d, maxncv, -6,
     &            'Ritz values and relative residuals')
          end if
c
c         %-------------------------------------------%
c         | Print additional convergence information. |
c         %-------------------------------------------%
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
          print *, ' _SSIMP '
          print *, ' ====== '
          print *, ' '
          print *, ' Size of the matrix is ', n
          print *, ' The number of Ritz values requested is ', nev
          print *, ' The number of Arnoldi vectors generated',
     &             ' (NCV) is ', ncv
          print *, ' What portion of the spectrum: ', which
          print *, ' The number of converged Ritz values is ', 
     &               nconv 
          print *, ' The number of Implicit Arnoldi update',
     &             ' iterations taken is ', iparam(3)
          print *, ' The number of OP*x is ', iparam(9)
          print *, ' The convergence criterion is ', tol
          print *, ' '
c
      end if
c
c     %---------------------------%
c     | Done with program sssimp. |
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
      Real
     &                  v(nx*nx), w(nx*nx), one, h2
      parameter         ( one = 1.0E+0 ) 
c
      call tv(nx,v(1),w(1))
      call saxpy(nx, -one, v(nx+1), 1, w(1), 1)
c
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call saxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
         call saxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
c
      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call saxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
c
c     Scale the vector w by (1/h^2), where h is the mesh size
c
      n2 = nx*nx
      h2 = one / real((nx+1)*(nx+1))
      call sscal(n2, one/h2, w, 1)
      return
      end
c
c-------------------------------------------------------------------
      subroutine tv (nx, x, y)
c
      integer           nx, j 
      Real
     &                  x(nx), y(nx), dd, dl, du
c
      Real
     &                  one, four
      parameter         (one = 1.0E+0, four = 4.0E+0)
c
c     Compute the matrix vector multiplication y<---T*x
c     where T is a nx by nx tridiagonal matrix with DD on the 
c     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
c     
c
      dd  = four 
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

