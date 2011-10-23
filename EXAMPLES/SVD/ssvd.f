      program ssvd
c
c     This example program is intended to illustrate the 
c     the use of ARPACK to compute the Singular Value Decomposition.
c   
c     This code shows how to use ARPACK to find a few of the
c     largest singular values(sigma) and corresponding right singular 
c     vectors (v) for the the matrix A by solving the symmetric problem:
c          
c                        (A'*A)*v = sigma*v
c 
c     where A is an m by n real matrix.
c
c     This code may be easily modified to estimate the 2-norm
c     condition number  largest(sigma)/smallest(sigma) by setting
c     which = 'BE' below.  This will ask for a few of the smallest
c     and a few of the largest singular values simultaneously.
c     The condition number could then be estimated by taking
c     the ratio of the largest and smallest singular values.
c
c     This formulation is appropriate when  m  .ge.  n.
c     Reverse the roles of A and A' in the case that  m .le. n.
c
c     The main points illustrated here are 
c
c        1) How to declare sufficient memory to find NEV 
c           largest singular values of A .  
c
c        2) Illustration of the reverse communication interface 
c           needed to utilize the top level ARPACK routine SSAUPD 
c           that computes the quantities needed to construct
c           the desired singular values and vectors(if requested).
c
c        3) How to extract the desired singular values and vectors
c           using the ARPACK routine SSEUPD.
c
c        4) How to construct the left singular vectors U from the 
c           right singular vectors V to obtain the decomposition
c
c                        A*V = U*S
c
c           where S = diag(sigma_1, sigma_2, ..., sigma_k).
c
c     The only thing that must be supplied in order to use this
c     routine on your problem is to change the array dimensions 
c     appropriately, to specify WHICH singular values you want to 
c     compute and to supply a the matrix-vector products 
c
c                         w <-  Ax
c                         y <-  A'w
c
c     in place of the calls  to AV( ) and ATV( ) respectively below.  
c
c     Further documentation is available in the header of DSAUPD
c     which may be found in the SRC directory.
c
c     This codes implements
c
c\Example-1
c     ... Suppose we want to solve A'A*v = sigma*v in regular mode,
c         where A is derived from the simplest finite difference 
c         discretization of the 2-dimensional kernel  K(s,t)dt  where
c
c                 K(s,t) =  s(t-1)   if 0 .le. s .le. t .le. 1,
c                           t(s-1)   if 0 .le. t .lt. s .le. 1. 
c
c         See subroutines AV  and ATV for details.
c     ... OP = A'*A  and  B = I.
c     ... Assume "call av (n,x,y)" computes y = A*x
c     ... Assume "call atv (n,y,w)" computes w = A'*y
c     ... Assume exact shifts are used
c     ...
c
c\BeginLib
c
c\Routines called:
c     ssaupd  ARPACK reverse communication interface routine.
c     sseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     sscal   Level 1 BLAS thst computes x <- x*alpha.
c     scopy   Level 1 BLAS thst computes y <- x.
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
c FILE: svd.F   SID: 2.4   DATE OF SID: 10/17/00   RELEASE: 2
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
c     | It is assumed that A is M by N with M .ge. N.        |
c     |                                                      |
c     | The maximum dimensions for all arrays are            |
c     | set here to accommodate a problem size of            |
c     | M .le. MAXM  and  N .le. MAXN                        |
c     |                                                      |
c     | The NEV right singular vectors will be computed in   |
c     | the N by NCV array V.                                |
c     |                                                      |
c     | The NEV left singular vectors will be computed in    |
c     | the M by NEV array U.                                |
c     |                                                      |
c     | NEV is the number of singular values requested.      |
c     |     See specifications for ARPACK usage below.       |
c     |                                                      |
c     | NCV is the largest number of basis vectors that will |
c     |     be used in the Implicitly Restarted Arnoldi      |
c     |     Process.  Work per major iteration is            |
c     |     proportional to N*NCV*NCV.                       |
c     |                                                      |
c     | You must set:                                        |
c     |                                                      |
c     | MAXM:   Maximum number of rows of the A allowed.     |
c     | MAXN:   Maximum number of columns of the A allowed.  |
c     | MAXNEV: Maximum NEV allowed                          |
c     | MAXNCV: Maximum NCV allowed                          |
c     %------------------------------------------------------%
c
      integer          maxm, maxn, maxnev, maxncv, ldv, ldu
      parameter       (maxm = 500, maxn=250, maxnev=10, maxncv=25, 
     &                 ldu = maxm, ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      Real
     &                 v(ldv,maxncv), u(ldu, maxnev), 
     &                 workl(maxncv*(maxncv+8)), workd(3*maxn), 
     &                 s(maxncv,2), resid(maxn), ax(maxm)
      logical          select(maxncv)
      integer          iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, m, n, nev, ncv, lworkl, info, ierr,
     &                 j, ishfts, maxitr, mode1, nconv
      logical          rvec
      Real      
     &                 tol, sigma, temp
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &                 one, zero
      parameter        (one = 1.0E+0, zero = 0.0E+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real           
     &                 snrm2
      external         snrm2, saxpy, scopy, sscal
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
      m = 500
      n = 100
c
c     %------------------------------------------------%
c     | Specifications for ARPACK usage are set        | 
c     | below:                                         |
c     |                                                |
c     |    1) NEV = 4 asks for 4 singular values to be |  
c     |       computed.                                | 
c     |                                                |
c     |    2) NCV = 20 sets the length of the Arnoldi  |
c     |       factorization                            |
c     |                                                |
c     |    3) This is a standard problem               |
c     |         (indicated by bmat  = 'I')             |
c     |                                                |
c     |    4) Ask for the NEV singular values of       |
c     |       largest magnitude                        |
c     |         (indicated by which = 'LM')            |
c     |       See documentation in SSAUPD for the      |
c     |       other options SM, BE.                    | 
c     |                                                |
c     | Note: NEV and NCV must satisfy the following   |
c     |       conditions:                              |
c     |                 NEV <= MAXNEV,                 |
c     |             NEV + 1 <= NCV <= MAXNCV           |
c     %------------------------------------------------%
c
      nev   = 4
      ncv   = 10 
      bmat  = 'I'
      which = 'LM'
c
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SVD: N is greater than MAXN '
         go to 9000
      else if ( m .gt. maxm ) then
         print *, ' ERROR with _SVD: M is greater than MAXM '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SVD: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SVD: NCV is greater than MAXNCV '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     | Specification of stopping rules and initial         |
c     | conditions before calling SSAUPD                    |
c     |                                                     |
c     |           abs(sigmaC - sigmaT) < TOL*abs(sigmaC)    |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |              (machine precision) is used.           |
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
c     | (indicated by setting IPARAM(1) = 1.)             |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of SSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | SSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = n
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
c           %---------------------------------------%
c           | Perform matrix vector multiplications |
c           |              w <--- A*x       (av())  |
c           |              y <--- A'*w      (atv()) |
c           | The user should supply his/her own    |
c           | matrix vector multiplication routines |
c           | here that takes workd(ipntr(1)) as    |
c           | the input, and returns the result in  |
c           | workd(ipntr(2)).                      |
c           %---------------------------------------%
c
            call av (m, n, workd(ipntr(1)), ax) 
            call atv (m, n, ax, workd(ipntr(2)))
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
c        %--------------------------------------------%
c        | No fatal errors occurred.                  |
c        | Post-Process using SSEUPD.                 |
c        |                                            |
c        | Computed singular values may be extracted. |  
c        |                                            |
c        | Singular vectors may also be computed now  |
c        | if desired.  (indicated by rvec = .true.)  | 
c        |                                            |
c        | The routine SSEUPD now called to do this   |
c        | post processing                            | 
c        %--------------------------------------------%
c           
         rvec = .true.
c
         call sseupd ( rvec, 'All', select, s, v, ldv, sigma, 
     &        bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
c
c        %-----------------------------------------------%
c        | Singular values are returned in the first     |
c        | column of the two dimensional array S         |
c        | and the corresponding right singular vectors  | 
c        | are returned in the first NEV columns of the  |
c        | two dimensional array V as requested here.    |
c        %-----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of SSEUPD. |
c           %------------------------------------%
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
               s(j,1) = sqrt(s(j,1))
c
c              %-----------------------------%
c              | Compute the left singular   |
c              | vectors from the formula    |
c              |                             |
c              |     u = Av/sigma            |
c              |                             |
c              | u should have norm 1 so     |
c              | divide by norm(Av) instead. |
c              %-----------------------------%
c
               call av(m, n, v(1,j), ax)
               call scopy(m, ax, 1, u(1,j), 1)
               temp = one/snrm2(m, u(1,j), 1)
               call sscal(m, temp, u(1,j), 1)
c
c              %---------------------------%
c              |                           |
c              | Compute the residual norm |
c              |                           |
c              |   ||  A*v - sigma*u ||    |
c              |                           |
c              | for the NCONV accurately  |
c              | computed singular values  |
c              | and vectors.  (iparam(5)  |
c              | indicates how many are    |
c              | accurate to the requested |
c              | tolerance).               |
c              | Store the result in 2nd   |
c              | column of array S.        |
c              %---------------------------%
c
               call saxpy(m, -s(j,1), u(1,j), 1, ax, 1)
               s(j,2) = snrm2(m, ax, 1)
c
 20         continue
c
c           %-------------------------------%
c           | Display computed residuals    |
c           %-------------------------------%
c
            call smout(6, nconv, 2, s, maxncv, -6,
     &                'Singular values and direct residuals')
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
         print *, ' _SVD '
         print *, ' ==== '
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
c     %-------------------------%
c     | Done with program ssvd. |
c     %-------------------------%
c
 9000 continue
c
      end
c 
c ------------------------------------------------------------------
c     matrix vector subroutines
c
c     The matrix A is derived from the simplest finite difference 
c     discretization of the integral operator 
c
c                     f(s) = integral(K(s,t)x(t)dt).
c      
c     Thus, the matrix A is a discretization of the 2-dimensional kernel 
c     K(s,t)dt, where
c
c                 K(s,t) =  s(t-1)   if 0 .le. s .le. t .le. 1,
c                           t(s-1)   if 0 .le. t .lt. s .le. 1.
c
c     Thus A is an m by n matrix with entries
c
c                 A(i,j) = k*(si)*(tj - 1)  if i .le. j,
c                          k*(tj)*(si - 1)  if i .gt. j
c
c     where si = i/(m+1)  and  tj = j/(n+1)  and k = 1/(n+1).
c      
c-------------------------------------------------------------------
c
      subroutine av (m, n, x, w)
c
c     computes  w <- A*x
c
      integer          m, n, i, j
      Real
     &                 x(n), w(m), one, zero, h, k, s, t
      parameter        ( one = 1.0E+0, zero = 0.0E+0 ) 
c
      h = one / real(m+1)
      k = one / real(n+1)
      do 5 i = 1,m
         w(i) = zero
 5    continue
      t = zero
c      
      do 30 j = 1,n
         t = t+k
         s = zero
         do 10 i = 1,j
           s = s+h
           w(i) = w(i) + k*s*(t-one)*x(j)
 10      continue 
         do 20 i = j+1,m
           s = s+h
           w(i) = w(i) + k*t*(s-one)*x(j) 
 20      continue
 30   continue      
c
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine atv (m, n, w, y)
c
c     computes  y <- A'*w
c
      integer         m, n, i, j
      Real
     &                w(m), y(n), one, zero,  h, k, s, t
      parameter       ( one = 1.0E+0, zero = 0.0E+0 )
c
      h = one / real(m+1)
      k = one / real(n+1)
      do 5 i = 1,n
         y(i) = zero
 5    continue
      t = zero
c
      do 30 j = 1,n
         t = t+k
         s = zero
         do 10 i = 1,j
           s = s+h
           y(j) = y(j) + k*s*(t-one)*w(i)
 10      continue
         do 20 i = j+1,m
           s = s+h
           y(j) = y(j) + k*t*(s-one)*w(i)
 20      continue
 30   continue
c
      return
      end 
c


