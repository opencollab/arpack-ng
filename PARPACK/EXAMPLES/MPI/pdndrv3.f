      program psndrv3
c
c     Message Passing Layer: MPI
c
c     Simple program to illustrate the idea of reverse communication
c     in inverse mode for a generalized nonsymmetric eigenvalue problem.
c
c     We implement example three of ex-nonsym.doc in DOCUMENTS directory
c
c\Example-3
c     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
c         where A is derived from the 1-dimensional convection-diffusion
c         operator on the interval [0,1] with zero boundary condition,
c         and M is the tridiagonal matrix with 4 on the diagonal and 1
c         on the subdiagonals.
c     ... So OP = inv[M]*A  and  B = M.
c     ... Use mode 2 of PDNAUPD.
c
c\BeginLib
c
c\Routines called:
c     pdnaupd Parallel ARPACK reverse communication interface routine.
c     pdneupd Parallel ARPACK routine that returns Ritz values and (optionally)
c              Ritz vectors.
c     dpttrf  LAPACK symmetric positive definite tridiagonal factorization
c             routine.
c     dpttrs  LAPACK symmetric positive definite tridiagonal solve routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     pdnorm2 Parallel version of Level 1 BLAS that computes the norm of a vector.
c     av      Parallel Matrix vector multiplication routine that computes A*x.
c     mv      Matrix vector multiplication routine that computes M*x.
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
c\Parallel Modifications
c     Kristi Maschhoff
c
c\Revision history:
c     Starting Point: Serial Code FILE: ndrv3.F   SID: 2.2
c
c\SCCS Information:
c FILE: ndrv3.F   SID: 1.1   DATE OF SID: 8/13/96   RELEASE: 1
c
c\Remarks
c     1. None
c
c\EndLib
c--------------------------------------------------------------------------
c
      include 'mpif.h'
      include 'debug.h'
      include 'stat.h'

c     %-------------------------------%
c     | MPI INTERFACE                 |
c     | ILP64 is not supported by MPI |
c     | integer*4 must be imposed in  |
c     | all calls involving MPI.      |
c     |                               |
c     | Use ierr for MPI calls.       |
c     %-------------------------------%

      integer*4         comm, myid, nprocs, rc, ierr

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
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Double precision
     &                  ax(maxn), mx(maxn), d(maxncv, 3), resid(maxn),
     &                  v(ldv,maxncv), workd(3*maxn),
     &                  workev(3*maxncv),
     &                  workl(3*maxncv*maxncv+6*maxncv),
     &                  md(maxn), me(maxn-1), temp(maxn), temp_buf(maxn)
c
c     %------------------------------------%
c     | Local Scalars                      |
c     |                                    |
c     | Use info if ILP64 can be supported |
c     | (call to BLAS, LAPACK, ARPACK).    |
c     %------------------------------------%
c
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, nloc, j,
     &                  nconv, maxitr, ishfts, mode, blk
      Double precision
     &                  tol, sigmar, sigmai
      logical           first, rvec
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                  zero, one
      parameter         (zero = 0.0, one = 1.0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
      Double precision
     &                  pdnorm2, dlapy2
      external          daxpy, pdnorm2, dpttrf, dpttrs, dlapy2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      call MPI_INIT( ierr )
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
c
      ndigit = -3
      logfil = 6
      mnaupd = 1
c
c     %----------------------------------------------------%
c     | The number N is the dimension of the matrix.  A    |
c     | generalized eigenvalue problem is solved (BMAT =   |
c     | 'G').  NEV is the number of eigenvalues to be      |
c     | approximated.  The user can modify NEV, NCV, WHICH |
c     | to solve problems of different sizes, and to get   |
c     | different parts of the spectrum.  However, The     |
c     | following conditions must be satisfied:            |
c     |                    N <= MAXN,                      |
c     |                  NEV <= MAXNEV,                    |
c     |              NEV + 2 <= NCV <= MAXNCV              |
c     %----------------------------------------------------%
c
      n     = 100
      nev   = 4
      ncv   = 20
c
c     %--------------------------------------%
c     | Set up distribution of data to nodes |
c     %--------------------------------------%
c
      nloc = (n / nprocs )
      blk = nloc
      if ( mod(n, nprocs) .gt. 0 ) then
         if ( myid .eq. nprocs-1 ) nloc = nloc + mod(n, nprocs)
*      if ( mod(n, nprocs) .gt. myid ) nloc = nloc + 1
      endif
c
      if ( nloc .gt. maxn ) then
         print *, ' ERROR with _NDRV3: NLOC is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
c
c     %-----------------------------------------------------%
c     | The matrix M is chosen to be the symmetric tri-     |
c     | diagonal matrix with 4 on the diagonal and 1 on the |
c     | off diagonals. It is factored by LAPACK subroutine  |
c     | dpttrf.                                             |
c     %-----------------------------------------------------%
c
      do 20 j = 1, n-1
         md(j) = 4.0
         me(j) = one
  20  continue
      md(n) = 4.0*one
c
      call dpttrf(n, md, me, info)
      if ( info .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _pttrf. '
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
      tol    = 0.0
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 2 of DNAUPD is used     |
c     | (IPARAM(7) = 2).  All these options can be        |
c     | changed by the user. For details, see the         |
c     | documentation in DNAUPD.                          |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 2
c
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
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
         call pdnaupd( comm, ido, bmat, nloc, which, nev, tol, resid,
     &                 ncv, v, ldv, iparam, ipntr, workd,
     &                 workl, lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %----------------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x      |
c           | The user should supply his/her own     |
c           | matrix vector routine and a linear     |
c           | system solver.  The matrix-vector      |
c           | subroutine should take workd(ipntr(1)) |
c           | as input, and the final result should  |
c           | be returned to workd(ipntr(2)).        |
c           %----------------------------------------%
c
            call av (comm, nloc, n, workd(ipntr(1)), workd(ipntr(2)))
c======== Hack for Linear system ======= ccc
            call dscal(n, zero, temp, 1)
            call dscal(n, zero, temp_buf, 1)
            do 15 j=1,nloc
               temp_buf(myid*blk + j) = workd(ipntr(2) + j - 1)
   15       continue
            call MPI_ALLREDUCE( temp_buf, temp, n,
     &            MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
            call dpttrs(n, 1, md, me, temp, n,
     &                  info)
            if ( info .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _pttrs. '
               print*, ' '
               go to 9000
            end if
            do 16 j=1,nloc
               workd(ipntr(2) + j - 1 ) = temp(myid*blk + j)
   16       continue
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         else if ( ido .eq. 2) then
c
c           %-------------------------------------%
c           |        Perform  y <--- M*x          |
c           | The matrix vector multiplication    |
c           | routine should take workd(ipntr(1)) |
c           | as input and return the result to   |
c           | workd(ipntr(2)).                    |
c           %-------------------------------------%
c
            call mv (comm, nloc, workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if
c
c
c     %-----------------------------------------%
c     | Either we have convergence, or there is |
c     | an error.                               |
c     %-----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %---------------------------%
c        | Error message. Check the  |
c        | documentation in PDNAUPD. |
c        %---------------------------%
c
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd.'
            print *, ' '
         endif
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
         call pdneupd ( comm, rvec, 'A', select, d, d(1,2), v, ldv,
     &        sigmar, sigmai, workev, bmat, nloc, which, nev, tol,
     &        resid, ncv, v, ldv, iparam, ipntr, workd,
     &        workl, lworkl, info )
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
         if ( info .ne. 0 ) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of DNEUPD. |
c           %------------------------------------%
c
            if ( myid .eq. 0 ) then
               print *, ' '
               print *, ' Error with _neupd, info = ', info
               print *, ' Check the documentation of _neupd'
               print *, ' '
            endif
c
         else
c
            first = .true.
            nconv = iparam(5)
            do 30 j=1, iparam(5)
c
c              %---------------------------%
c              | Compute the residual norm |
c              |                           |
c              |  ||  A*x - lambda*M*x ||  |
c              |                           |
c              | for the NCONV accurately  |
c              | computed eigenvalues and  |
c              | eigenvectors.  (iparam(5) |
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
                  call av(comm, nloc, n, v(1,j), ax)
                  call mv(comm, nloc, v(1,j), mx)
                  call daxpy(nloc, -d(j,1), mx, 1, ax, 1)
                  d(j,3) = pdnorm2(comm, nloc, ax, 1)
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
                  call av(comm, nloc, n, v(1,j), ax)
                  call mv(comm, nloc, v(1,j), mx)
                  call daxpy(nloc, -d(j,1), mx, 1, ax, 1)
                  call mv(comm, nloc, v(1,j+1), mx)
                  call daxpy(nloc, d(j,2), mx, 1, ax, 1)
                  d(j,3) = pdnorm2(comm, nloc, ax, 1)**2
                  call av(comm, nloc, n, v(1,j+1), ax)
                  call mv(comm, nloc, v(1,j+1), mx)
                  call daxpy(nloc, -d(j,1), mx, 1, ax, 1)
                  call mv(comm, nloc, v(1,j), mx)
                  call daxpy(nloc, -d(j,2), mx, 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), pdnorm2(comm,nloc,ax,1) )
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
c
  30        continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call pdmout(comm, 6, nconv, 3, d, maxncv, -6,
     &           'Ritz values (Real,Imag) and direct residuals')
c
         end if
c
c        %------------------------------------------%
c        | Print additional convergence information |
c        %------------------------------------------%
c
         if (myid .eq. 0)then
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit
     &                 Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
c
         print *, ' '
         print *, '_NDRV3 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processors is ', nprocs
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
         endif
      end if
c
c     %----------------------------%
c     | Done with program psndrv3. |
c     %----------------------------%
c
 9000 continue
c
c     %-------------------------%
c     | Release resources MPI |
c     %-------------------------%
c
      call MPI_FINALIZE(rc)
c
      end
c
c==========================================================================
c
c     parallel matrix vector multiplication subroutine
c
c     Compute the matrix vector multiplication y<---A*x
c     where A is a n by n nonsymmetric tridiagonal matrix derived
c     from the central difference discretization of the 1-dimensional
c     convection diffusion operator on the interval [0,1] with
c     zero Dirichlet boundary condition.
c
      subroutine av (comm, nloc, n, v, w)
c
c     .. MPI Declarations ...
      include           'mpif.h'
      integer*4         comm, nprocs, myid, ierr,
     &                  status(MPI_STATUS_SIZE)
c
      integer           nloc, n, j, next, prev
      Double precision
     &                  v(nloc), w(nloc), one, two, dd, dl, du,
     &                  s, h, rho, mv_buf
      parameter         ( rho = 10.0, one = 1.0,
     &                    two = 2.0)
c
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
      h = one / dble(n+1)
      s = rho*h / two
      dd = two
      dl = -one - s
      du = -one + s
c
      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,nloc-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(nloc) =  dl*v(nloc-1) + dd*v(nloc)
c
      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v(nloc), 1, MPI_DOUBLE_PRECISION,
     &                  next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, prev, myid,
     &                  comm, status, ierr )
         w(1) = w(1) + dl*mv_buf
      endif
c
      if ( myid .gt. 0 ) then
         call mpi_send( v(1), 1, MPI_DOUBLE_PRECISION,
     &                  prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, next, myid,
     &                  comm, status, ierr )
         w(nloc) = w(nloc) + du*mv_buf
      endif
c
      return
      end
c------------------------------------------------------------------------
c
c     Compute the matrix vector multiplication y<---M*x
c     where M is a n by n tridiagonal matrix with 4 on the
c     diagonal, 1 on the subdiagonal and the superdiagonal.
c
      subroutine mv (comm, nloc, v, w)
c
c     .. MPI Declarations ...
      include           'mpif.h'
      integer*4         comm, nprocs, myid, ierr,
     &                  status(MPI_STATUS_SIZE)
c
      integer           nloc, j, next, prev
      Double precision
     &                  v(nloc), w(nloc), one, four, mv_buf
      parameter         ( one = 1.0, four = 4.0)
c
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
c
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,nloc-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1)
 10   continue
      w(nloc) =  one*v(nloc-1) + four*v(nloc)
c
      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v(nloc), 1, MPI_DOUBLE_PRECISION,
     &                  next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, prev, myid,
     &                  comm, status, ierr )
         w(1) = w(1) + mv_buf
      endif
c
      if ( myid .gt. 0 ) then
         call mpi_send( v(1), 1, MPI_DOUBLE_PRECISION,
     &                  prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, next, myid,
     &                  comm, status, ierr )
         w(nloc) = w(nloc) + mv_buf
      endif
c
      return
      end
c------------------------------------------------------------
      subroutine mv2 (comm, n, v, w)
      integer           n, j, comm
      Double precision
     &                  v(n), w(n)
      do 10 j=1,n
         w(j) = v(j)
 10   continue
c
      return
      end
