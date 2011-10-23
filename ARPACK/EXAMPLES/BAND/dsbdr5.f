      program dsbdr5 
c
c     ... Construct the matrix A in LAPACK-style band form.
c         The matrix A is the 1-dimensional discrete Laplacian on [0,1]
c         with zero Dirichlet boundary condition, KG is the mass
c         formed by using piecewise linear elements on [0,1]. 
c
c     ... Call DSBAND  with Buckling mode to find eigenvalues LAMBDA 
c         such that
c                          A*x = M*x*LAMBDA.
c
c     ... Use mode 4 of DSAUPD .
c
c\BeginLib
c
c\Routines called:
c     dsband   ARPACK banded eigenproblem solver.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlaset   LAPACK routine to initialize a matrix to zero.
c     daxpy    Level 1 BLAS that computes y <- alpha*x+y.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dgbmv    Level 2 BLAS that computes the band matrix vector product
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
c FILE: sbdr5.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2
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
      Double precision 
     &                 a(lda,maxn), m(lda,maxn), rfac(lda,maxn),
     &                 workl(maxncv*maxncv+8*maxncv), workd(3*maxn), 
     &                 v(ldv, maxncv), resid(maxn), d(maxncv, 2),
     &                 ax(maxn), mx(maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        which*2, bmat
      integer          nev, ncv, kl, ku, info, j, ido,
     &                 n, isub, isup, idiag, maxitr, mode, nconv
      Double precision  
     &                 tol, h, sigma, r1, r2
      logical          rvec
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision  
     &                 one, zero, two, four, six
      parameter        (one = 1.0D+0 , zero = 0.0D+0 , two = 2.0D+0 ,
     &                  four = 4.0D+0 , six = 6.0D+0 )
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision 
     &                  dlapy2 , dnrm2 
      external          dlapy2 , dnrm2 , dgbmv , daxpy  
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
c     %--------------------------------------------------%
c     | The number N is the dimension of the matrix.  A  |
c     | generalized eigenvalue problem is solved         |
c     | (BMAT = 'G').  NEV is the number of eigenvalues  |
c     | to be approximated.  Since the Buckling mode is  |
c     | is used, WHICH is set to 'LM'.  The user can     |
c     | modify N, NEV, NCV and SIGMA to solve problems   |
c     | of different sizes, and to get different parts   |
c     | the spectrum.  However, the following conditions |
c     | must be satisfied:                               |
c     |                   N <= MAXN                      |
c     |                 NEV <= MAXNEV                    |
c     |           NEV + 1 <= NCV <= MAXNCV               | 
c     %--------------------------------------------------% 
c
      n    = 100
      nev  = 4 
      ncv  = 10 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _SBDR5: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _SBDR5: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _SBDR5: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'G'
      which = 'LM'
      sigma = 1.0
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in DSAUPD  as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in DSAUPD  to start the Arnoldi iteration. |
c     %-----------------------------------------------------%
c
      lworkl  = ncv**2+8*ncv
      tol  = zero 
      ido  = 0
      info = 0
c
c     %---------------------------------------------------%
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 4 of DSAUPD  is used     |
c     | (IPARAM(7) = 4). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DSBAND .                                           |
c     %---------------------------------------------------%
c
      maxitr = 300
      mode   = 4
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
      call dlaset ('A', lda, n, zero, zero, a, lda)
      call dlaset ('A', lda, n, zero, zero, m, lda)
      call dlaset ('A', lda, n, zero, zero, rfac, lda)
c
c     %-------------------------------------%
c     | KU, KL are number of superdiagonals |
c     | and subdiagonals within the band of |
c     | matrices A and M.                   |
c     %-------------------------------------%
c
      kl   = 1 
      ku   = 1 
c
c     %---------------% 
c     | Main diagonal |
c     %---------------%
c
      h = one / dble (n+1)
      r1 = four / six
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = two / h
         m(idiag,j) = r1 * h
  30  continue 
c 
c     %-------------------------------------%
c     | First subdiagonal and superdiagonal |
c     %-------------------------------------%
c 
      r2 = one / six
      isup = kl+ku
      isub = kl+ku+2
      do 60 j = 1, n-1
         a(isup,j+1) = -one / h
         a(isub,j) = -one / h
         m(isup,j+1) =  r2 * h
         m(isub,j) = r2 * h
  60  continue
c
c     %-------------------------------------%
c     | Call DSBAND  to find eigenvalues and |
c     | eigenvectors.  Eigenvalues are      |
c     | returned in the first column of D.  |
c     | Eigenvectors are returned in the    |
c     | first NCONV (=IPARAM(5)) columns of |
c     | V.                                  |
c     %-------------------------------------%
c
      rvec = .true.
      call dsband ( rvec, 'A', select, d, v, ldv, sigma, n, a, m, lda, 
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
         print *, ' _SBDR5 '
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
            call dgbmv ('Notranspose', n, n, kl, ku, one, 
     &                 a(kl+1,1), lda, v(1,j), 1, zero, 
     &                 ax, 1)
            call dgbmv ('Notranspose', n, n, kl, ku, one, 
     &                 m(kl+1,1), lda, v(1,j), 1, zero, 
     &                 mx, 1)
            call daxpy (n, -d(j,1), mx, 1, ax, 1)
            d(j,2) = dnrm2 (n, ax, 1)
            d(j,2) = d(j,2) / abs(d(j,1))
c
 90      continue 

         call dmout (6, nconv, 2, d, maxncv, -6,
     &             'Ritz values and relative residuals')
      else 
c
c        %-------------------------------------%
c        | Either convergence failed, or there |
c        | is error.  Check the documentation  |
c        | for DSBAND .                         |
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
