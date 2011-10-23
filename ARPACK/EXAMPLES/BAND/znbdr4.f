      program zndrv4 
c
c     ... Construct matrices A and M in LAPACK-style band form.
c         Matries A and M are derived from the finite
c         element discretization of the 1-dimensional
c         convection-diffusion operator
c                         (d^2u/dx^2) + rho*(du/dx)
c         on the interval [0,1] with zero boundary condition using
c         piecewise linear elements.
c
c     ... Call ZNBAND  to find eigenvalues LAMBDA such that
c                    A*x = M*x*LAMBDA.
c
c     ... Use mode 3 of ZNAUPD .
c
c\BeginLib
c
c\Routines called:
c     znband   ARPACK banded eigenproblem solver.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     zlaset   LAPACK routine to initialize a matrix to zero.
c     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
c     dznrm2   Level 1 BLAS that computes the norm of a vector.
c     zgbmv    Level 2 BLAS that computes the band matrix vector product.
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
c FILE: nbdr4.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-------------------------------------------------------------------------
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
     &                   maxbdw=50, lda = maxbdw, ldv = maxn)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          iparam(11), iwork(maxn)
      logical          select(maxncv)
      Complex*16  
     &                 a(lda,maxn), m(lda,maxn), fac(lda,maxn),
     &                 workl(3*maxncv*maxncv+5*maxncv), workd(3*maxn), 
     &                 workev(2*maxncv), v(ldv, maxncv),
     &                 resid(maxn), d(maxncv), ax(maxn), mx(maxn)
      Double precision  
     &                 rwork(maxn), rd(maxncv,3)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, j,
     &                 n, idiag, isup, isub, maxitr, mode,
     &                 nconv
      logical          rvec
      Double precision  
     &                 tol
      Complex*16 
     &                 rho, h, sigma
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16  
     &                 one, zero, two, four, six
      parameter        (one  = (1.0D+0, 0.0D+0)  ,
     &                  zero = (0.0D+0, 0.0D+0) ,
     &                  two  = (2.0D+0, 0.0D+0) ,
     &                  four = (4.0D+0, 0.0D+0) ,
     &                  six  = (6.0D+0, 0.0D+0) )
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision 
     &                  dznrm2 , dlapy2 
      external          dznrm2 , zgbmv , zaxpy , dlapy2 , zlaset 
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------------------------%
c     | The number NX is the number of interior points  |
c     | in the discretization of the 2-dimensional      |
c     | convection-diffusion operator on the unit       |
c     | square with zero Dirichlet boundary condition.  |
c     | The number N(=NX*NX) is the dimension of the    |
c     | matrix.  A generalized eigenvalue problem is    |
c     | solved (BMAT = 'G').  NEV is the number of      |
c     | eigenvalues (closest to the shift SIGMA) to be  |
c     | approximated.  Since the shift and invert mode  |
c     | is used, WHICH is set to 'LM'.  The user can    |
c     | modify NX, NEV and NCV to solve problems of     | 
c     | different sizes, and to get different parts the |
c     | spectrum.  However, the following conditions    |
c     | must be satisfied:                              |
c     |                   N <= MAXN                     |
c     |                 NEV <= MAXNEV                   |
c     |           NEV + 2 <= NCV <= MAXNCV              |
c     %-------------------------------------------------%
c
      n    = 100
      nev  = 4 
      ncv  = 10 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NBDR4: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NBDR4: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NBDR4: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'G'
      which = 'LM'
      sigma = (1.0D+1, 0.0D+0) 
c
c     %----------------------------------------------------%
c     | The work array WORKL is used in ZNAUPD  as          | 
c     | workspace.  Its dimension LWORKL has to be set as  |
c     | illustrated below.  The parameter TOL determines   |
c     | the stopping criterion. If TOL<=0, machine machine |
c     | precision is used.  Setting INFO=0 indicates that  |
c     | we using a randomly generated vector to start the  |
c     | the ARNOLDI process.                               | 
c     %----------------------------------------------------%
c
      lworkl  = 3*ncv**2+5*ncv
      info = 0
      tol  = 0.0
c
c     %---------------------------------------------------%
c     | IPARAm(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 3 of ZNAUPD  is used     |
c     | (IPARAm(7) = 3). All these options can be changed |
c     | by the user. For details, see the documentation   |
c     | in znband .                                        |
c     %---------------------------------------------------%
c
      maxitr = 300
      mode   = 3
c
      iparam(3) = maxitr
      iparam(7) = mode
c
c     %--------------------------------------------%
c     | Construct matrices A and M in LAPACK-style |
c     | banded form.                               |
c     %--------------------------------------------%
c
c     %---------------------------------------------%
c     | Zero out the workspace for banded matrices. |
c     %---------------------------------------------%
c
      call zlaset ('A', lda, n, zero, zero, a, lda)
      call zlaset ('A', lda, n, zero, zero, m, lda)
      call zlaset ('A', lda, n, zero, zero, fac, lda)
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
      h = one / dcmplx (n+1)
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = two / h 
         m(idiag,j) = four * h / six
  30  continue 
c 
c     %-------------------------------------%
c     | First subdiagonal and superdiagonal |
c     %-------------------------------------%
c 
      isup = kl+ku
      isub = kl+ku+2
      rho = (1.0D+1, 0.0D+0) 
      do 40 j = 1, n-1
         a(isup,j+1) = -one/h + rho/two
         a(isub,j) = -one/h - rho/two
         m(isup,j+1) = one*h / six
         m(isub,j) = one*h / six
  40  continue 
c
c     %-----------------------------------------------%
c     | Call ARPACK banded solver to find eigenvalues |
c     | and eigenvectors. Eigenvalues are returned in |
c     | the one dimensional array D.  Eigenvectors    |
c     | are returned in the first NCONV (=IPARAM(5))  |
c     | columns of V.                                 |
c     %-----------------------------------------------%
c
      rvec = .true. 
      call znband (rvec, 'A', select, d, v, ldv, sigma,
     &           workev, n, a, m, lda, fac, kl, ku, which,
     &           bmat, nev, tol, resid, ncv, v, ldv, iparam,
     &           workd, workl, lworkl, rwork, iwork, info)
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
         print *, '_NBDR4 '
         print *, '====== '
         print *, ' '
         print *, ' The size of the matrix is ', n
         print *, ' Number of eigenvalue requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     &            ' (NCV) is ', ncv
         print *, ' The number of converged Ritz values is ',
     &              nconv
         print *, ' What portion of the spectrum ', which
         print *, ' The number of Implicit Arnoldi ',
     &              ' update taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence tolerance is ', tol
         print *, ' '
c
         do 90 j = 1, nconv
c
c           %----------------------------%
c           | Compute the residual norm. |
c           |    ||  A*x - lambda*x ||   |
c           %----------------------------%
c
            call zgbmv ('Notranspose', n, n, kl, ku, one,
     &                 a(kl+1,1), lda, v(1,j), 1, zero,
     &                 ax, 1)
            call zgbmv ('Notranspose', n, n, kl, ku, one,
     &                 m(kl+1,1), lda, v(1,j), 1, zero,
     &                 mx, 1)
            call zaxpy (n, -d(j), mx, 1, ax, 1)
            rd(j,1) = dble (d(j))
            rd(j,2) = dimag (d(j))
            rd(j,3) = dznrm2 (n, ax, 1)
            rd(j,3) = rd(j,3) / dlapy2 (rd(j,1), rd(j,2))
 90      continue 

         call dmout (6, nconv, 3, rd, maxncv, -6,
     &             'Ritz values (Real,Imag) and relative residuals')
      else 
c
c        %-------------------------------------%
c        | Either convergence failed, or there |
c        | is error.  Check the documentation  |
c        | for znband .                         |
c        %-------------------------------------%
c
          print *, ' '
          print *, ' Error with _band, info= ', info
          print *, ' Check the documentation of _band '
          print *, ' ' 
c
      end if
c
 9000 end      
