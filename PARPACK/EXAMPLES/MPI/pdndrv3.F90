      program psndrv3
!
!     Message Passing Layer: MPI
!
!     Simple program to illustrate the idea of reverse communication
!     in inverse mode for a generalized nonsymmetric eigenvalue problem.
!
!     We implement example three of ex-nonsym.doc in DOCUMENTS directory
!
!\Example-3
!     ... Suppose we want to solve A*x = lambda*B*x in inverse mode,
!         where A is derived from the 1-dimensional convection-diffusion
!         operator on the interval [0,1] with zero boundary condition,
!         and M is the tridiagonal matrix with 4 on the diagonal and 1
!         on the subdiagonals.
!     ... So OP = inv[M]*A  and  B = M.
!     ... Use mode 2 of PDNAUPD.
!
!\BeginLib
!
!\Routines called:
!     pdnaupd Parallel ARPACK reverse communication interface routine.
!     pdneupd Parallel ARPACK routine that returns Ritz values and (optionally)
!              Ritz vectors.
!     dpttrf  LAPACK symmetric positive definite tridiagonal factorization
!             routine.
!     dpttrs  LAPACK symmetric positive definite tridiagonal solve routine.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
!     pdnorm2 Parallel version of Level 1 BLAS that computes the norm of a vector.
!     av      Parallel Matrix vector multiplication routine that computes A*x.
!     mv      Matrix vector multiplication routine that computes M*x.
!
!\Author
!     Richard Lehoucq
!     Danny Sorensen
!     Chao Yang
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Code FILE: ndrv3.F   SID: 2.2
!
!\SCCS Information:
! FILE: ndrv3.F   SID: 1.1   DATE OF SID: 8/13/96   RELEASE: 1
!
!\Remarks
!     1. None
!
!\EndLib
!--------------------------------------------------------------------------
!
#ifdef HAVE_MPI_ICB
      use :: mpi_f08
#else
#include "mpif.h"
#endif
#include "debugF90.h"
#include "statF90.h"

!     %-------------------------------%
!     | MPI INTERFACE                 |
!     | ILP64 is not supported by MPI |
!     | integer*4 must be imposed in  |
!     | all calls involving MPI.      |
!     | MPI communicators must be     |
!     | declared with MPI_Comm type.  |
!     | Use mpi_f08 to get correct    |
!     | types (= ICB provided by MPI) |
!     |                               |
!     | Use ierr for MPI calls.       |
!     %-------------------------------%

#ifdef HAVE_MPI_ICB
      type(MPI_Comm)    comm
#else
      integer*4         comm
#endif
      integer*4         myid, nprocs, rc, ierr

!
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=10, maxncv=25, ldv=maxn )
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Double precision  ax(maxn), mx(maxn), d(maxncv, 3), resid(maxn),&
                        v(ldv,maxncv), workd(3*maxn),&
                        workev(3*maxncv),&
                        workl(3*maxncv*maxncv+6*maxncv),&
                        md(maxn), me(maxn-1), temp(maxn), temp_buf(maxn)
!
!     %------------------------------------%
!     | Local Scalars                      |
!     |                                    |
!     | Use info if ILP64 can be supported |
!     | (call to BLAS, LAPACK, ARPACK).    |
!     %------------------------------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j,&
                        nconv, maxitr, ishfts, mode, blk
      integer*4         nloc
      Double precision  tol, sigmar, sigmai
      logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Double precision  zero, one
      parameter         (zero = 0.0, one = 1.0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
      Double precision  pdnorm2, dlapy2
      external          daxpy, pdnorm2, dpttrf, dpttrs, dlapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      call MPI_INIT( ierr )
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
!
      ndigit = -3
      logfil = 6
      mnaupd = 1
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues to be      |
!     | approximated.  The user can modify NEV, NCV, WHICH |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum.  However, The     |
!     | following conditions must be satisfied:            |
!     |                    N <= MAXN,                      |
!     |                  NEV <= MAXNEV,                    |
!     |              NEV + 2 <= NCV <= MAXNCV              |
!     %----------------------------------------------------%
!
      n     = 100
      nev   = 4
      ncv   = 20
!
!     %--------------------------------------%
!     | Set up distribution of data to nodes |
!     %--------------------------------------%
!
      nloc = (n / nprocs )
      blk = nloc
      if ( mod(n, nprocs) .gt. 0 ) then
         if ( myid .eq. nprocs-1 ) nloc = nloc + mod(n, nprocs)
      endif
!
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
!
!     %-----------------------------------------------------%
!     | The matrix M is chosen to be the symmetric tri-     |
!     | diagonal matrix with 4 on the diagonal and 1 on the |
!     | off diagonals. It is factored by LAPACK subroutine  |
!     | dpttrf.                                             |
!     %-----------------------------------------------------%
!
      do 20 j = 1, n-1
         md(j) = 4.0
         me(j) = one
  20  continue
      md(n) = 4.0*one
!
      call dpttrf(n, md, me, info)
      if ( info .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _pttrf. '
         print*, ' '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in DNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in DNAUPD to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
!
      lworkl = 3*ncv**2+6*ncv
      tol    = 0.0
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 2 of DNAUPD is used     |
!     | (IPARAM(7) = 2).  All these options can be        |
!     | changed by the user. For details, see the         |
!     | documentation in DNAUPD.                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 2
!
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call pdnaupd( comm, ido, bmat, nloc, which, nev, tol, resid,&
                       ncv, v, ldv, iparam, ipntr, workd,&
                       workl, lworkl, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %----------------------------------------%
!           | Perform  y <--- OP*x = inv[M]*A*x      |
!           | The user should supply his/her own     |
!           | matrix vector routine and a linear     |
!           | system solver.  The matrix-vector      |
!           | subroutine should take workd(ipntr(1)) |
!           | as input, and the final result should  |
!           | be returned to workd(ipntr(2)).        |
!           %----------------------------------------%
!
            call av (comm, nloc, n, workd(ipntr(1)), workd(ipntr(2)))
!======== Hack for Linear system ======= ccc
            call dscal(n, zero, temp, 1)
            call dscal(n, zero, temp_buf, 1)
            do 15 j=1,nloc
               temp_buf(myid*blk + j) = workd(ipntr(2) + j - 1)
   15       continue
            call MPI_ALLREDUCE( temp_buf, temp, n,&
                  MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
            call dpttrs(n, 1, md, me, temp, n,&
                        info)
            if ( info .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _pttrs. '
               print*, ' '
               go to 9000
            end if
            do 16 j=1,nloc
               workd(ipntr(2) + j - 1 ) = temp(myid*blk + j)
   16       continue
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         else if ( ido .eq. 2) then
!
!           %-------------------------------------%
!           |        Perform  y <--- M*x          |
!           | The matrix vector multiplication    |
!           | routine should take workd(ipntr(1)) |
!           | as input and return the result to   |
!           | workd(ipntr(2)).                    |
!           %-------------------------------------%
!
            call mv (comm, nloc, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         end if
!
!
!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %---------------------------%
!        | Error message. Check the  |
!        | documentation in PDNAUPD. |
!        %---------------------------%
!
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd.'
            print *, ' '
         endif
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
         call pdneupd ( comm, rvec, 'A', select, d, d(1,2), v, ldv,&
              sigmar, sigmai, workev, bmat, nloc, which, nev, tol,&
              resid, ncv, v, ldv, iparam, ipntr, workd,&
              workl, lworkl, info )
!
!        %-----------------------------------------------%
!        | The real part of the eigenvalue is returned   |
!        | in the first column of the two dimensional    |
!        | array D, and the IMAGINARY part is returned   |
!        | in the second column of D.  The corresponding |
!        | eigenvectors are returned in the first NEV    |
!        | columns of the two dimensional array V if     |
!        | requested.  Otherwise, an orthogonal basis    |
!        | for the invariant subspace corresponding to   |
!        | the eigenvalues in D is returned in V.        |
!        %-----------------------------------------------%
!
         if ( info .ne. 0 ) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
!
            if ( myid .eq. 0 ) then
               print *, ' '
               print *, ' Error with _neupd, info = ', info
               print *, ' Check the documentation of _neupd'
               print *, ' '
            endif
!
         else
!
            first = .true.
            nconv = iparam(5)
            do 30 j=1, iparam(5)
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |  ||  A*x - lambda*M*x ||  |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (iparam(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               if (d(j,2) .eq. zero)  then
!
!                 %--------------------%
!                 | Ritz value is real |
!                 %--------------------%
!
                  call av(comm, nloc, n, v(1,j), ax)
                  call mv(comm, nloc, v(1,j), mx)
                  call daxpy(nloc, -d(j,1), mx, 1, ax, 1)
                  d(j,3) = pdnorm2(comm, nloc, ax, 1)
!
               else if (first) then
!
!                 %------------------------%
!                 | Ritz value is complex  |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
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
!
  30        continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call pdmout(comm, 6, nconv, 3, d, maxncv, -6,&
                 'Ritz values (Real,Imag) and direct residuals')
!
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if (myid .eq. 0)then
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' '
            print *, ' No shifts could be applied during implicit'&
                   , ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if
!
         print *, ' '
         print *, '_NDRV3 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processors is ', nprocs
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',&
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',&
                    nconv
         print *, ' The number of Implicit Arnoldi update',&
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
         endif
      end if
!
!     %----------------------------%
!     | Done with program psndrv3. |
!     %----------------------------%
!
 9000 continue
!
!     %-------------------------%
!     | Release resources MPI |
!     %-------------------------%
!
      call MPI_FINALIZE(rc)
!
      end
!
!==========================================================================
!
!     parallel matrix vector multiplication subroutine
!
!     Compute the matrix vector multiplication y<---A*x
!     where A is a n by n nonsymmetric tridiagonal matrix derived
!     from the central difference discretization of the 1-dimensional
!     convection diffusion operator on the interval [0,1] with
!     zero Dirichlet boundary condition.
!
      subroutine av (comm, nloc, n, v, w)
!
!     .. MPI Declarations ...
#ifdef HAVE_MPI_ICB
      use :: mpi_f08
#else
#include "mpif.h"
#endif
#ifdef HAVE_MPI_ICB
      type(MPI_Comm)    comm
      type(MPI_Status)  status
#else
      integer*4         comm, status(MPI_STATUS_SIZE)
#endif
      integer*4         nprocs, myid, ierr, next, prev
!
      integer*4         nloc, n, j
      Double precision  v(nloc), w(nloc), one, two, dd, dl, du,&
                        s, h, rho, mv_buf
      parameter         ( rho = 10.0, one = 1.0, two = 2.0)
!
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
      h = one / dble(n+1)
      s = rho*h / two
      dd = two
      dl = -one - s
      du = -one + s
!
      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,nloc-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(nloc) =  dl*v(nloc-1) + dd*v(nloc)
!
      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v(nloc), 1, MPI_DOUBLE_PRECISION,&
                        next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, prev, myid,&
                        comm, status, ierr )
         w(1) = w(1) + dl*mv_buf
      endif
!
      if ( myid .gt. 0 ) then
         call mpi_send( v(1), 1, MPI_DOUBLE_PRECISION,&
                        prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, next, myid,&
                        comm, status, ierr )
         w(nloc) = w(nloc) + du*mv_buf
      endif
!
      return
      end
!------------------------------------------------------------------------
!
!     Compute the matrix vector multiplication y<---M*x
!     where M is a n by n tridiagonal matrix with 4 on the
!     diagonal, 1 on the subdiagonal and the superdiagonal.
!
      subroutine mv (comm, nloc, v, w)
!
!     .. MPI Declarations ...
#ifdef HAVE_MPI_ICB
      use :: mpi_f08
#else
#include "mpif.h"
#endif
#ifdef HAVE_MPI_ICB
      type(MPI_Comm)    comm
      type(MPI_Status)  status
#else
      integer*4         comm, status(MPI_STATUS_SIZE)
#endif
      integer*4         nprocs, myid, ierr, next, prev
!
      integer*4         nloc, j
      Double precision  v(nloc), w(nloc), one, four, mv_buf
      parameter         ( one = 1.0, four = 4.0)
!
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
!
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,nloc-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1)
 10   continue
      w(nloc) =  one*v(nloc-1) + four*v(nloc)
!
      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v(nloc), 1, MPI_DOUBLE_PRECISION,&
                        next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, prev, myid,&
                        comm, status, ierr )
         w(1) = w(1) + mv_buf
      endif
!
      if ( myid .gt. 0 ) then
         call mpi_send( v(1), 1, MPI_DOUBLE_PRECISION,&
                        prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, 1, MPI_DOUBLE_PRECISION, next, myid,&
                        comm, status, ierr )
         w(nloc) = w(nloc) + mv_buf
      endif
!
      return
      end
!------------------------------------------------------------
      subroutine mv2 (comm, n, v, w)
      integer           n, j, comm
      Double precision  v(n), w(n)
      do 10 j=1,n
         w(j) = v(j)
 10   continue
!
      return
      end
