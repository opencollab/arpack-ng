      program snbdr3
c
c     ... Construct matrices A and M in LAPACK-style band form.
c         The matrix A and M are derived from the finite element
c         discretization of the 1-dimensional convection-diffusion operator 
c                         (d^2u/dx^2) + rho*(du/dx)
c         on the interval [0,1] with zero boundary condition,

c     ... Call SNBAND to find eigenvalues LAMBDA such that
c                    A*x = LAMBDA*M*x.
c
c     ... Eigenvalues with largest real parts are sought.
c
c     ... Use mode 2 of SNAUPD.
c
c\BeginLib
c
c\Routines called:
c     snband  ARPACK banded eigenproblem solver.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     slaset  LAPACK routine to initialize a matrix to zero.
c     saxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     sgbmv   Level 2 BLAS that computes the band matrix vector product.
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
c FILE: nbdr3.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2
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
      Real 
     &                 a(lda,maxn), m(lda,maxn), rfac(lda,maxn),
     &                 workl(3*maxncv*maxncv+6*maxncv), workd(3*maxn), 
     &                 workev(3*maxncv), v(ldv, maxncv),
     &                 resid(maxn), d(maxncv, 3), ax(maxn), mx(maxn)
      Complex  
     &                 cfac(lda, maxn), workc(maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, j, ido,
     &                 n, idiag, isup, isub, mode, maxitr,
     &                 nconv
      logical          rvec, first
      Real  
     &                 tol, rho, h, sigmar, sigmai
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Real  
     &                 one, zero, two
      parameter        (one = 1.0E+0 , zero = 0.0E+0 , 
     &                  two = 2.0E+0 )
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
c     | The number N is the dimension of the matrix.  A |
c     | generalized eigenvalue problem is solved        |
c     | (BMAT = 'G').  NEV is the number of eigenvalues |
c     | to be approximated. The user can modify N, NEV, |
c     | NCV and WHICH to solve problems of different    |
c     | sizes, and to get different parts the spectrum. |
c     | However, the following conditions must be       |
c     | satisfied:                                      |
c     |                   N <= MAXN                     |
c     |                 NEV <= MAXNEV                   |
c     |           NEV + 2 <= NCV <= MAXNCV              |
c     %-------------------------------------------------%
c
      n    = 100
      nev  = 4 
      ncv  = 10 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NBDR3: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NBDR3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NBDR3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'G'
      which = 'LM'
c
c     %----------------------------------------------------%
c     | The work array WORKL is used in SNAUPD as          | 
c     | workspace.  Its dimension LWORKL has to be set as  |
c     | illustrated below.  The parameter TOL determines   |
c     | the stopping criterion. If TOL<=0, machine machine |
c     | precision is used.  The number IDO is used for     |
c     | reverse communication and has to be set to 0 at    |
c     | the beginning.  Setting INFO=0 indicates that we   |
c     | using a randomly generated vector to start the     |
c     | the ARNOLDI process.                               | 
c     %----------------------------------------------------%
c
      lworkl  = 3*ncv**2+6*ncv
      info = 0
      tol  = zero 
      ido  = 0
c
c     %---------------------------------------------------%
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 2 of SNAUPD is used     |
c     | (IPARAM(7) = 2). All these options can be changed |
c     | by the user. For details, see the documentation   |
c     | in SNBAND.                                        |
c     %---------------------------------------------------%
c
      mode = 2
      maxitr = 300
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
      kl   = 1 
      ku   = 1 
c
c     %---------------% 
c     | Main diagonal |
c     %---------------%
c
      h = one / real (n+1)
c
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = 2.0E+0  / h
         m(idiag,j) = 4.0E+0  * h
  30  continue 
c 
c     %-------------------------------------%
c     | First subdiagonal and superdiagonal |
c     %-------------------------------------%
c 
      isup = kl+ku
      isub = kl+ku+2
      rho = 1.0E+1 
      do 50 j = 1, n
         a(isup,j+1) = -one/h + rho/two
         a(isub,j) = -one/h - rho/two
         m(isup,j+1) = one*h
         m(isub,j) = one*h
  50  continue 
c
c     %------------------------------------------------%
c     | Call ARPACK banded solver to find eigenvalues  |
c     | and eigenvectors. The real parts of the        |
c     | eigenvalues are returned in the first column   |
c     | of D, the imaginary parts are returned in the  |
c     | second column of D.  Eigenvectors are returned |
c     | in the first NCONV (=IPARAM(5)) columns of V.  |
c     %------------------------------------------------%
c
      rvec = .true. 
      call snband( rvec, 'A', select, d, d(1,2), v, ldv, sigmar, 
     &     sigmai, workev, n, A, M, lda, rfac, cfac, kl, ku, 
     &     which, bmat, nev, tol, resid, ncv, v, ldv, iparam, 
     &     workd, workl, lworkl, workc, iwork, info)
c
      if ( info .eq. 0) then
c
c        %-----------------------------------%
c        | Print out convergence information |
c        %-----------------------------------%
c
         nconv = iparam(5)
c
         print *, ' '
         print *, ' _NBDR3 '
         print *, ' ====== '
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
c        %----------------------------%
c        | Compute the residual norm. |
c        |    ||  A*x - lambda*x ||   |
c        %----------------------------%
c
         first = .true. 
         do 90 j = 1, nconv
c
            if ( d(j,2) .eq. zero ) then
c
c              %--------------------%
c              | Ritz value is real |
c              %--------------------%
c
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    a(kl+1,1), lda, v(1,j), 1, zero, 
     &                    ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    m(kl+1,1), lda, v(1,j), 1, zero, 
     &                    mx, 1)
               call saxpy(n, -d(j,1), mx, 1, ax, 1)
               d(j,3) = snrm2(n, ax, 1)
               d(j,3) = d(j,3) / abs(d(j,1))
c
            else if ( first ) then
c
c              %------------------------%
c              | Ritz value is complex  |
c              | Residual of one Ritz   |
c              | value of the conjugate |
c              | pair is computed.      | 
c              %------------------------%
c
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    a(kl+1,1), lda, v(1,j), 1, zero, 
     &                    ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    m(kl+1,1), lda, v(1,j), 1, zero, 
     &                    mx, 1)
               call saxpy(n, -d(j,1), mx, 1, ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    m(kl+1,1), lda, v(1,j+1), 1, zero, 
     &                    mx, 1)
               call saxpy(n, d(j,2), mx, 1, ax, 1)
               d(j,3) = snrm2(n, ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    a(kl+1,1), lda, v(1,j+1), 1, zero, 
     &                    ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    m(kl+1,1), lda, v(1,j+1), 1, zero, 
     &                    mx, 1)
               call saxpy(n, -d(j,1), mx, 1, ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    m(kl+1,1), lda, v(1,j), 1, zero, 
     &                    mx, 1)
               call saxpy(n, -d(j,2), mx, 1, ax, 1)
               d(j,3) = slapy2( d(j,3), snrm2(n, ax, 1) )
               d(j,3) = d(j,3) / slapy2(d(j,1),d(j,2))
               d(j+1,3) = d(j,3)
               first = .false.
            else
               first = .true.
            end if
c
 90      continue 

         call smout(6, nconv, 3, d, maxncv, -6,
     &             'Ritz values (Real,Imag) and relative residuals')
      else 
c
c        %-------------------------------------%
c        | Either convergence failed, or there |
c        | is error.  Check the documentation  |
c        | for SNBAND.                         |
c        %-------------------------------------%
c
          print *, ' '
          print *, ' Error with _nband, info= ', info
          print *, ' Check the documentation of _nband '
          print *, ' ' 
c
      end if
c
 9000 end      
