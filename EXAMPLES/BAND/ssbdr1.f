      program ssbdr1
c
c     ... Construct the matrix A in LAPACK-style band form.
c         The matrix A is derived from the discretization of
c         the 2-dimensional Laplacian on the unit square with 
c         zero Dirichlet boundary condition using standard 
c         central difference.
c
c     ... Call SSBAND to find eigenvalues LAMBDA such that
c                          A*x = x*LAMBDA.
c       
c     ... Use mode 1 of SSAUPD.
c
c\BeginLib
c
c\Routines called:
c     ssband  ARPACK banded eigenproblem solver.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     slaset  LAPACK routine to initialize a matrix to zero.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     sgbmv   Level 2 BLAS that computes the band matrix vector product
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
c FILE: sbdr1.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c----------------------------------------------------------------------
c
c     %-------------------------------------%
c     | Define leading dimensions for all   |
c     | arrays.                             |
c     | MAXN   - Maximum size of the matrix |
c     | MAXNEV - Maximum number of          |
c     |          eigenvalues to be computed |
c     | MAXNCV - Maximum number of Arnoldi  |
c     |          vectors stored             | 
c     | MAXBDW - Maximum bandwidth          |
c     %-------------------------------------%
c
      integer          maxn, maxnev, maxncv, maxbdw, lda,
     &                 lworkl, ldv
      parameter        ( maxn = 1000, maxnev = 25, maxncv=50, 
     &                   maxbdw=50, lda = maxbdw, ldv = maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          iparam(11), iwork(maxn)
      logical          select(maxncv)
      Real 
     &                 a(lda,maxn), m(lda,maxn), rfac(lda,maxn),
     &                 workl(maxncv*maxncv+8*maxncv), workd(3*maxn), 
     &                 v(ldv, maxncv), resid(maxn), d(maxncv, 2),
     &                 ax(maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, i, j, ido,
     &                 n, nx, lo, isub, isup, idiag, maxitr, mode,
     &                 nconv
      Real  
     &                 tol, sigma, h2
      logical          rvec
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Real  
     &                 one, zero, two
      parameter        (one = 1.0E+0 , zero = 0.0E+0 , two = 2.0E+0 )
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real 
     &                  slapy2, snrm2
      external          slapy2, snrm2, sgbmv, saxpy 
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
c     | The number NX is the number of interior points  |
c     | in the discretization of the 2-dimensional      |
c     | Laplacian operator on the unit square with zero |
c     | Dirichlet boundary condition. The number        |
c     | N(=NX*NX) is the dimension of the matrix.  A    |
c     | standard eigenvalue problem is solved           |
c     | (BMAT = 'I').  NEV is the number of eigenvalues |
c     | to be approximated. The user can modify NX,NEV, |
c     | NCV and WHICH to solve problems of different    |
c     | sizes, and to get different parts the spectrum. |
c     | However, the following conditions must be       |
c     | satisfied:                                      |
c     |                   N <= MAXN                     |
c     |                 NEV <= MAXNEV                   |
c     |           NEV + 1 <= NCV <= MAXNCV              | 
c     %-------------------------------------------------% 
c
      nx  = 10 
      n    = nx*nx
      nev  = 4 
      ncv  = 10 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SBDR1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SBDR1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SBDR1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'LM'
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in SSAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in SSAUPD to start the Arnoldi iteration. |
c     %-----------------------------------------------------%
c
      lworkl  = ncv**2+8*ncv
      tol  = zero 
      ido  = 0
      info = 0
c
c     %---------------------------------------------------%
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of SSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | SSBAND.                                           |
c     %---------------------------------------------------%
c
      maxitr  = 300
      mode = 1
c
      iparam(3) = maxitr
      iparam(7) = mode 
c
c     %----------------------------------------%
c     | Construct the matrix A in LAPACK-style |
c     | banded form.                           |
c     %----------------------------------------%
c
c     %---------------------------------------------%
c     | Zero out the workspace for banded matrices. |
c     %---------------------------------------------%
c
      call slaset('A', lda, n, zero, zero, a, lda)
      call slaset('A', lda, n, zero, zero, m, lda)
      call slaset('A', lda, n, zero, zero, rfac, lda)
c
c     %-------------------------------------%
c     | KU, KL are number of superdiagonals |
c     | and subdiagonals within the band of |
c     | matrices A and M.                   |
c     %-------------------------------------%
c
      kl   = nx 
      ku   = nx 
c
c     %---------------% 
c     | Main diagonal |
c     %---------------%
c
      h2 = one / ((nx+1)*(nx+1))
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = 4.0E+0  / h2
  30  continue 
c 
c     %-------------------------------------%
c     | First subdiagonal and superdiagonal |
c     %-------------------------------------%
c 
      isup = kl+ku
      isub = kl+ku+2
      do 50 i = 1, nx
        lo = (i-1)*nx
        do 40 j = lo+1, lo+nx-1
           a(isup,j+1) = -one / h2
           a(isub,j) = -one / h2
  40    continue      
  50  continue 
c
c     %------------------------------------%
c     | KL-th subdiagonal and KU-th super- |
c     | diagonal.                          |
c     %------------------------------------%
c
      isup = kl+1
      isub = 2*kl+ku+1
      do 80 i = 1, nx-1
         lo = (i-1)*nx
         do 70 j = lo+1, lo+nx
            a(isup,nx+j)  = -one / h2
            a(isub,j) = -one / h2
 70      continue 
 80   continue 
c
c     %-------------------------------------%
c     | Call SSBAND to find eigenvalues and |
c     | eigenvectors.  Eigenvalues are      |
c     | returned in the first column of D.  |
c     | Eigenvectors are returned in the    |
c     | first NCONV (=IPARAM(5)) columns of |
c     | V.                                  |
c     %-------------------------------------% 
c
      rvec = .true.
      call ssband( rvec, 'A', select, d, v, ldv, sigma, n, a, m, lda, 
     &             rfac, kl, ku, which, bmat, nev, tol, 
     &             resid, ncv, v, ldv, iparam, workd, workl, lworkl, 
     &             iwork, info)
c
      if ( info .eq. 0) then
c
         nconv = iparam(5)
c 
c        %-----------------------------------%
c        | Print out convergence information |
c        %-----------------------------------%
c
         print *, ' '
         print *, ' _SBDR1 '
         print *, ' ====== '
         print *, ' '
         print *, ' The size of the matrix is ', n
         print *, ' Number of eigenvalue requested is ', nev
         print *, ' The number of Lanczos vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' The number of converged Ritz values is ',
     &              nconv
         print *, ' What portion of the spectrum ', which
         print *, ' The number of Implicit Arnoldi',
     &              ' update taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence tolerance is ', tol
         print *, ' '
c
c        %----------------------------%
c        | Compute the residual norm. |
c        |    ||  A*x - lambda*x ||   |
c        %----------------------------%
c
         do 90 j = 1, nconv
            call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                 a(kl+1,1), lda, v(1,j), 1, zero, 
     &                 ax, 1)
            call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
            d(j,2) = snrm2(n, ax, 1)
            d(j,2) = d(j,2) / abs(d(j,1))
c
 90      continue 

         call smout(6, nconv, 2, d, maxncv, -6,
     &             'Ritz values and relative residuals')
      else 
c
c        %-------------------------------------%
c        | Either convergence failed, or there |
c        | is error.  Check the documentation  |
c        | for SSBAND.                         |
c        %-------------------------------------%
c
          print *, ' '
          print *, ' Error with _sband, info= ', info
          print *, ' Check the documentation of _sband '
          print *, ' ' 
c
      end if
c
 9000 end      
