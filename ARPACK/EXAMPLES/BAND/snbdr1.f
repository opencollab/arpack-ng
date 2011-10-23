      program snbdr1
c
c     ... Construct the matrix A in LAPACK-style band form.
c         The matrix A is derived from the discretization of
c         the 2-d convection-diffusion operator
c
c              -Laplacian(u) + rho*partial(u)/partial(x).
c
c         on the unit square with zero Dirichlet boundary condition
c         using standard central difference.
c
c     ... Call SNBAND to find eigenvalues LAMBDA such that
c                          A*x = LAMBDA*x.
c
c     ... Use mode 1 of SNAUPD.
c
c\BeginLib
c
c     snband  ARPACK banded eigenproblem solver.
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
c FILE: nbdr1.F   SID: 2.5   DATE OF SID: 08/26/96   RELEASE: 2
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
     &                 workl(3*maxncv*maxncv+6*maxncv), workd(3*maxn), 
     &                 workev(3*maxncv), v(ldv, maxncv),
     &                 resid(maxn), d(maxncv, 3), ax(maxn)
      Complex  
     &                 cfac(lda, maxn), workc(maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        which*2, bmat
      integer          nev, ncv, ku, kl, info, i, j, ido,
     &                 n, nx, lo, isub, isup, idiag, mode, maxitr,
     &                 nconv
      logical          rvec, first
      Real  
     &                 tol, rho, h, h2, sigmar, sigmai
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
c     | The number NX is the number of interior points  |
c     | in the discretization of the 2-dimensional      |
c     | convection-diffusion operator on the unit       |
c     | square with zero Dirichlet boundary condition.  | 
c     | The number N(=NX*NX) is the dimension of the    |
c     | matrix.  A standard eigenvalue problem is       |
c     | solved (BMAT = 'I').  NEV is the number of      |
c     | eigenvalues to be approximated. The user can    |
c     | modify NX, NEV, NCV, WHICH to solve problems of |
c     | different sizes, and to get different parts the |
c     | spectrum.  However, The following conditions    |
c     | must be satisfied:                              |
c     |                   N <= MAXN                     |
c     |                 NEV <= MAXNEV                   |
c     |           NEV + 2 <= NCV <= MAXNCV              | 
c     %-------------------------------------------------% 
c
      nx  = 10 
      n    = nx*nx
      nev  = 4 
      ncv  = 10 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NBDR1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NBDR1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NBDR1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat = 'I'
      which = 'SM'
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in SNAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in SNAUPD to start the Arnoldi iteration. |
c     %-----------------------------------------------------%
c
      lworkl  = 3*ncv**2+6*ncv
      tol  = zero 
      ido  = 0
      info = 0
c
c     %---------------------------------------------------%
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of SNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details, see the documentation   |
c     | in SNBAND.                                        |
c     %---------------------------------------------------%
c
      maxitr = 300
      mode   = 1
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
      h  = one / real (nx+1)
      h2 = h*h
c
      idiag = kl+ku+1
      do 30 j = 1, n
         a(idiag,j) = 4.0E+0  / h2
  30  continue 
c 
c     %-------------------------------------%
c     | First subdiagonal and superdiagonal |
c     %-------------------------------------%
c 
      rho = 1.0E+2  
      isup = kl+ku
      isub = kl+ku+2
      do 50 i = 1, nx
        lo = (i-1)*nx
        do 40 j = lo+1, lo+nx-1
           a(isup,j+1) = -one/h2 + rho/two/h
           a(isub,j) = -one/h2 - rho/two/h
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
      call snband(rvec, 'A', select, d, d(1,2), v, ldv, sigmar, sigmai, 
     &           workev, n, a, m, lda, rfac, cfac, kl, ku, which, 
     &           bmat, nev, tol, resid, ncv, v, ldv, iparam, workd, 
     &           workl, lworkl, workc, iwork, info)
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
         print *, ' _NBDR1 '
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
               call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
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
               call saxpy(n, -d(j,1), v(1,j), 1, ax, 1)
               call saxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
               d(j,3) = snrm2(n, ax, 1)
               call sgbmv('Notranspose', n, n, kl, ku, one, 
     &                    a(kl+1,1), lda, v(1,j+1), 1, zero, 
     &                    ax, 1)
               call saxpy(n, -d(j,2), v(1,j), 1, ax, 1)
               call saxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
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
