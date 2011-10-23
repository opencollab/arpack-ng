      program dndrv6
c
c     Simple program to illustrate the idea of reverse communication
c     in shift-invert mode for a generalized nonsymmetric eigenvalue problem.
c
c     We implement example six of ex-nonsym.doc in DOCUMENTS directory
c
c\Example-6
c
c     ... Suppose we want to solve A*x = lambda*B*x in shift-invert mode
c         The matrix A is the tridiagonal matrix with 2 on the diagonal, 
c         -2 on the subdiagonal and 3 on the superdiagonal.  The matrix M 
c         is the tridiagonal matrix with 4 on the diagonal and 1 on the 
c         off-diagonals.
c     ... The shift sigma is a complex number (sigmar, sigmai).
c     ... OP = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M  and  B = M.
c     ... Use mode 4 of DNAUPD.
c
c\BeginLib
c
c\Routines called:
c     dnaupd  ARPACK reverse communication interface routine.
c     dneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     zgttrf  LAPACK complex matrix factorization routine.
c     zgttrs  LAPACK complex linear system solve routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     av      Matrix vector subroutine that computes A*x.
c     mv      Matrix vector subroutine that computes M*x.
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
c FILE: ndrv6.F   SID: 2.5   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c-------------------------------------------------------------------------
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
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=10, maxncv=25, 
     &                   ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14), ipiv(maxn)
      logical           select(maxncv)
      Double precision
     &                  ax(maxn), mx(maxn), d(maxncv,3), resid(maxn),
     &                  v(ldv,maxncv), workd(3*maxn),
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv)
      Complex*16         
     &                  cdd(maxn), cdl(maxn), cdu(maxn),
     &                  cdu2(maxn), ctemp(maxn)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, ierr, j,
     &                  nconv, maxitr, ishfts, mode
      Double precision
     &                  tol, numr, numi, denr, deni, sigmar, sigmai
      Complex*16         
     &                  c1, c2, c3
      logical           first, rvec 
c 
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      external          zgttrf, zgttrs
      Double precision   
     &                  ddot, dnrm2, dlapy2
      external          ddot, dnrm2, dlapy2
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                   zero
      parameter          (zero = 0.0D+0)
c
c     %--------------------%
c     | Intrinsic Function |
c     %--------------------%
c
      intrinsic          dimag, dcmplx, abs
c
c     %-----------------------%
c     | Executable statements |
c     %-----------------------%
c
c     %----------------------------------------------------%
c     | The number N is the dimension of the matrix.  A    |
c     | generalized eigenvalue problem is solved (BMAT =   |
c     | 'G').  NEV is the number of eigenvalues (closest   |
c     | to the shift (SIGMAR,SIGMAI)) to be approximated.  |
c     | Since the shift-invert mode is used, WHICH is set  |
c     | to 'LM'.  The user can modify NEV, NCV, SIGMA to   |
c     | solve problems of different sizes, and to get      |
c     | different parts of the spectrum. However, The      |
c     | following conditions must be satisfied:            |
c     |                     N <= MAXN,                     | 
c     |                   NEV <= MAXNEV,                   |
c     |               NEV + 2 <= NCV <= MAXNCV             | 
c     %----------------------------------------------------%
c
      n     = 100 
      nev   = 4 
      ncv   = 20 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV6: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV6: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV6: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
      sigmar = 4.0D-1
      sigmai = 6.0D-1 
c
c     %----------------------------------------------------%
c     | Construct C = A - (SIGMAR,SIGMAI)*M in complex     |
c     | arithmetic, and factor C in complex arithmetic     |
c     | (using LAPACK subroutine zgttrf). The matrix A is  |
c     | chosen to be the tridiagonal matrix with -2 on the |
c     | subdiagonal, 2 on the diagonal and 3 on the        |
c     | superdiagonal. The matrix M is chosen to be the    |
c     | symmetric tridiagonal matrix with 4 on the         |
c     | diagonal and 1 on the off-diagonals.               |
c     %----------------------------------------------------%
c
      c1 = dcmplx(-2.0D+0-sigmar, -sigmai)
      c2 = dcmplx( 2.0D+0-4.0D+0*sigmar, -4.0D+0*sigmai)
      c3 = dcmplx( 3.0D+0-sigmar, -sigmai)
c
      do 10 j = 1, n-1
         cdl(j) = c1 
         cdd(j) = c2
         cdu(j) = c3
  10  continue 
      cdd(n) = c2 
c 
      call zgttrf(n, cdl, cdd, cdu, cdu2, ipiv, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf in _NDRV6.'
         print*, ' '
         go to 9000
      end if
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in DNAUPD as           |
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in DNAUPD to start the Arnoldi iteration. |
c     %-----------------------------------------------------%
c
      lworkl = 3*ncv**2+6*ncv 
      tol    = zero 
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shift with respect to     |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 3 of DNAUPD is used     |
c     | (IPARAM(7) = 3).  All these options can be        |
c     | changed by the user. For details, see the         |
c     | documentation in DNAUPD.                          |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 3
c
      iparam(1) = ishfts
      iparam(3) = maxitr  
      iparam(7) = mode
c
c     %------------------------------------------%
c     | M A I N   L O O P(Reverse communication) | 
c     %------------------------------------------%
c
 20   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DNAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dnaupd ( ido, bmat, n, which, nev, tol, resid, 
     &                 ncv, v, ldv, iparam, ipntr, workd, 
     &                 workl, lworkl, info )
c

         if (ido .eq. -1) then
c
c           %------------------------------------------------------------%
c           |                           Perform                          | 
c           | y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} |
c           | to force starting vector into the range of OP. The user    |
c           | should supply his/her own matrix vector multiplication     |
c           | routine and a complex linear system solver.  The matrix    |
c           | vector multiplication routine should take workd(ipntr(1))  |
c           | as the input. The final result (a real vector) should be   |
c           | returned to workd(ipntr(2)).                               |
c           %------------------------------------------------------------%
c
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
            do 30 j = 1, n
               ctemp(j) = dcmplx(workd(ipntr(2)+j-1))
  30        continue
c
            call zgttrs('N', n, 1, cdl, cdd, cdu, cdu2, ipiv, 
     &                  ctemp, maxn, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV6.'
               print*, ' '
               go to 9000
            end if 
            do  40 j = 1, n
               workd(ipntr(2)+j-1) = dimag(ctemp(j))
  40        continue
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 20
c
         else if ( ido .eq. 1) then
c
c           %------------------------------------------------------------%
c           |                          Perform                           |
c           | y <--- OP*x = Imaginary_Part{inv[A-(SIGMAR,SIGMAI)*M]*M*x} |
c           | M*x has been saved in workd(ipntr(3)). The user only need  |
c           | the complex linear system solver here that takes           |
c           | complex[workd(ipntr(3))] as input, and returns the result  |
c           | to workd(ipntr(2)).                                        |
c           %------------------------------------------------------------%
c
            do 50 j = 1,n
               ctemp(j) = dcmplx(workd(ipntr(3)+j-1))
  50        continue
            call zgttrs ('N', n, 1, cdl, cdd, cdu, cdu2, ipiv, 
     &                   ctemp, maxn, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV6.'
               print*, ' '
               go to 9000
            end if 
            do  60 j = 1, n
               workd(ipntr(2)+j-1) = dimag(ctemp(j))
  60        continue
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 20
c
         else if ( ido .eq. 2) then
c
c           %---------------------------------------------%
c           |          Perform  y <--- M*x                |
c           | Need matrix vector multiplication routine   |
c           | here that takes workd(ipntr(1)) as input    |
c           | and returns the result to workd(ipntr(2)).  |
c           %---------------------------------------------%
c
            call mv (n, workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 20
c
         end if
c
c     %-----------------------------------------%
c     | Either we have convergence, or there is |
c     | an error.                               |
c     %-----------------------------------------%
c
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
         print *, ' Check the documentation of _naupd.'
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
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        %-------------------------------------------%
c
         rvec = .true.
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv,  
     &        sigmar, sigmai, workev, bmat, n, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, 
     &        workl, lworkl, ierr )
c
c        %-----------------------------------------------%
c        | The real part of the eigenvalue is returned   |
c        | in the first column of the two dimensional    |
c        | array D, and the IMAGINARY part is returned   |
c        | in the second column of D.  The corresponding |
c        | eigenvectors are returned in the first NEV    |
c        | columns of the two dimensional array V if     |
c        | requested.  Otherwise, an orthogonal basis    |
c        | for the invariant subspace corresponding to   |
c        | the eigenvalues in D is returned in V.        |
c        %-----------------------------------------------%
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
            do 70 j=1, nconv
c
c               %-------------------------------------%
c               | Use Rayleigh Quotient to recover    |
c               | eigenvalues of the original problem.|
c               %-------------------------------------%
c
                if ( d(j,2) .eq. zero) then
c
c                  %----------------------------%
c                  | Eigenvalue is real.        | 
c                  | Compute d = x'(Ax)/x'(Mx). |
c                  %----------------------------%
c
                   call av(n, v(1,j), ax )
                   numr = ddot(n, v(1,j), 1, ax, 1)
                   call mv(n, v(1,j), ax )
                   denr = ddot(n, v(1,j), 1, ax, 1)
                   d(j,1) =  numr / denr 
c
                else if (first) then
c
c                  %------------------------%
c                  | Eigenvalue is complex. |
c                  | Compute the first one  |
c                  | of the conjugate pair. |
c                  %------------------------%
c
c                  %----------------%
c                  | Compute x'(Ax) |
c                  %----------------%
c
                   call av(n, v(1,j), ax )
                   numr = ddot(n, v(1,j), 1, ax, 1)
                   numi = ddot(n, v(1,j+1), 1, ax, 1)
                   call av(n, v(1,j+1), ax)
                   numr = numr + ddot(n,v(1,j+1),1,ax,1) 
                   numi = -numi + ddot(n,v(1,j),1,ax,1)
c
c                  %----------------%
c                  | Compute x'(Mx) |
c                  %----------------%
c
                   call mv(n, v(1,j), ax )
                   denr = ddot(n, v(1,j), 1, ax, 1)
                   deni = ddot(n, v(1,j+1), 1, ax, 1)
                   call mv(n, v(1,j+1), ax)
                   denr = denr + ddot(n,v(1,j+1),1,ax,1) 
                   deni = -deni + ddot(n,v(1,j),1, ax,1)
c
c                  %----------------%
c                  | d=x'(Ax)/x'(Mx)|
c                  %----------------%
c
                   d(j,1) = (numr*denr+numi*deni) /
     &                      dlapy2(denr, deni)
                   d(j,2) = (numi*denr-numr*deni) /
     &                      dlapy2(denr, deni)
                   first = .false.
c
                else
c
c                  %------------------------------%
c                  | Get the second eigenvalue of |
c                  | the conjugate pair by taking |
c                  | the conjugate of the last    |
c                  | eigenvalue computed.         |
c                  %------------------------------%
c 
                   d(j,1) = d(j-1,1)
                   d(j,2) = -d(j-1,2)
                   first = .true.
c
                end if
c
  70        continue
c
c           %---------------------------%
c           | Compute the residual norm |
c           |                           |
c           |   ||  A*x - lambda*x ||   |
c           |                           |
c           | for the NCONV accurately  |
c           | computed eigenvalues and  |
c           | eigenvectors.  (iparam(5) |
c           | indicates how many are    |
c           | accurate to the requested |
c           | tolerance)                |
c           %---------------------------%
c
            first = .true.
            nconv = iparam(5)
            do 80 j=1, nconv 
c
               if (d(j,2) .eq. zero)  then
c
c                 %--------------------%
c                 | Ritz value is real |
c                 %--------------------%
c
                  call av(n, v(1,j), ax)
                  call mv(n, v(1,j), mx)
                  call daxpy(n, -d(j,1), mx, 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
c
               else if (first) then
c
c                 %------------------------%
c                 | Ritz value is complex  |
c                 | Residual of one Ritz   |
c                 | value of the conjugate |
c                 | pair is computed.      | 
c                 %------------------------%
c        
                  call av(n, v(1,j), ax)
                  call mv(n, v(1,j), mx)
                  call daxpy(n, -d(j,1), mx, 1, ax, 1)
                  call mv(n, v(1,j+1), mx)
                  call daxpy(n, d(j,2), mx, 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  call av(n, v(1,j+1), ax)
                  call mv(n, v(1,j+1), mx)
                  call daxpy(n, -d(j,1), mx, 1, ax, 1)
                  call mv(n, v(1,j), mx)
                  call daxpy(n, -d(j,2), mx, 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2)) 
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
c
  80        continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call dmout(6, nconv, 3, d, maxncv, -6,
     &           'Ritz values (Real,Imag) and relative residuals')
c
        end if
c
c       %-------------------------------------------%
c       | Print additional convergence information. |
c       %-------------------------------------------%
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
         print *, ' _NDRV6 '
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
c     | Done with program dndrv6. |
c     %---------------------------%
c
 9000 continue
c
      end
c 
c==========================================================================
c
c     matrix vector multiplication subroutine
c
      subroutine mv (n, v, w)
      integer           n, j
      Double precision
     &                  v(n), w(n), one, four 
      parameter         (one = 1.0D+0, four = 4.0D+0)
c
c     Compute the matrix vector multiplication y<---M*x
c     where M is a n by n symmetric tridiagonal matrix with 4 on the 
c     diagonal, 1 on the subdiagonal and superdiagonal.
c 
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,n-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1) 
 10   continue 
      w(n) =  one*v(n-1) + four*v(n) 
      return
      end
c------------------------------------------------------------------
      subroutine av (n, v, w)
      integer           n, j
      Double precision            
     &                  v(n), w(n), three, two 
      parameter         (three = 3.0D+0, two = 2.0D+0)
c
c     Compute the matrix vector multiplication y<---A*x
c     where M is a n by n symmetric tridiagonal matrix with 2 on the
c     diagonal, -2 on the subdiagonal and 3 on the superdiagonal.
c
      w(1) =  two*v(1) + three*v(2)
      do 10 j = 2,n-1
         w(j) = -two*v(j-1) + two*v(j) + three*v(j+1)
 10   continue
      w(n) =  -two*v(n-1) + two*v(n)
      return
      end

