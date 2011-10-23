      program psndrv1 
c
c     Message Passing Layer: MPI
c
c     Example program to illustrate the idea of reverse communication
c     for a standard nonsymmetric eigenvalue problem.
c
c     We implement example one of ex-nonsym.doc in DOCUMENTS directory
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is obtained from the standard central difference
c         discretization of the convection-diffusion operator 
c                 (Laplacian u) + rho*(du / dx)
c         on the unit square, with zero Dirichlet boundary condition.
c
c     ... OP = A  and  B = I.
c     ... Assume "call av (comm, nloc, nx, mv_buf, x, y)" computes y = A*x.
c     ... Use mode 1 of PSNAUPD.
c
c\BeginLib
c
c\Routines called:
c     psnaupd  Parallel ARPACK reverse communication interface routine.
c     psneupd  Parallel ARPACK routine that returns Ritz values and (optionally)
c              Ritz vectors.
c     slapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     saxpy    Level 1 BLAS that computes y <- alpha*x+y.
c     psnorm2  Parallel version of Level 1 BLAS that computes the norm of a vector.
c     av       Distributed matrix vector multiplication routine that computes A*x.
c     tv       Matrix vector multiplication routine that computes T*x, 
c              where T is a tridiagonal matrix.  It is used in routine
c              av.
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
c     Starting Point: Serial Code FILE: ndrv1.F   SID: 2.2
c
c\SCCS Information: 
c FILE: ndrv1.F   SID: 1.2   DATE OF SID: 8/14/96   RELEASE: 1
c
c\Remarks
c     1. None
c
c\EndLib
c---------------------------------------------------------------------------
c
      include 'mpif.h'
      include 'debug.h'
      include 'stat.h'
 
c     %---------------%
c     | MPI INTERFACE |
c     %---------------%
      integer           comm, myid, nprocs, rc, nloc 
c
c     %-----------------------------%
c     | Define maximum dimensions   |
c     | for all arrays.             |
c     | MAXN:   Maximum dimension   |
c     |         of the distributed  |
c     |         block of A allowed. |
c     | MAXNEV: Maximum NEV allowed |
c     | MAXNCV: Maximum NCV allowed |
c     %-----------------------------%
c
      integer           maxnloc, maxnev, maxncv, ldv
      parameter         (maxnloc=256, maxnev=12, maxncv=30, ldv=maxnloc)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Real
     &                  ax(maxnloc), d(maxncv,3), resid(maxnloc), 
     &                  v(ldv,maxncv), workd(3*maxnloc), 
     &                  workev(3*maxncv), 
     &                  workl(3*maxncv*maxncv+6*maxncv)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
      Real
     &                  tol, sigmar, sigmai
      logical           first, rvec
c
c     %----------------------------------------------%
c     | Local Buffers needed for MPI communication |
c     %----------------------------------------------%
c
      Real
     &                  mv_buf(maxnloc)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &                  zero
      parameter         (zero = 0.0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real
     &                  slapy2, psnorm2
      external          slapy2, saxpy, psnorm2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic         abs, sqrt
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
c     %--------------------------------------------------%
c     | The number NX is the number of interior points   |
c     | in the discretization of the 2-dimensional       |
c     | convection-diffusion operator on the unit        |
c     | square with zero Dirichlet boundary condition.   | 
c     | The number N(=NX*NX) is the dimension of the     |
c     | matrix.  A standard eigenvalue problem is        |
c     | solved (BMAT = 'I').  NEV is the number of       |
c     | eigenvalues to be approximated.  The user can    |
c     | modify NX, NEV, NCV, WHICH to solve problems of  |
c     | different sizes, and to get different parts of   |
c     | the spectrum.  However, The following            |
c     | conditions must be satisfied:                    |
c     |                   N <= MAXN                      |
c     |                 NEV <= MAXNEV                    |
c     |           NEV + 2 <= NCV <= MAXNCV               | 
c     %--------------------------------------------------% 
c
      nx    = 10 
      n     = nx*nx 
      nev   = 4
      ncv   = 20 
c
c     %--------------------------------------%
c     | Set up distribution of data to nodes |
c     %--------------------------------------%
c
      nloc = (nx / nprocs)*nx
      if ( mod(nx, nprocs) .gt. myid ) nloc = nloc + nx
c
      if ( nloc .gt. maxnloc ) then
         print *, ' ERROR with _NDRV1: NLOC is greater than MAXNLOC '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'SM'
c
c     %-----------------------------------------------------%
c     | The work array WORKL is used in PSNAUPD as          |  
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.  The parameter TOL determines    |
c     | the stopping criterion. If TOL<=0, machine          |
c     | precision is used.  The variable IDO is used for    |
c     | reverse communication, and is initially set to 0.   |
c     | Setting INFO=0 indicates that a random vector is    |
c     | generated in PSNAUPD to start the Arnoldi iteration.| 
c     %-----------------------------------------------------%
c
      lworkl  = 3*ncv**2+6*ncv 
      tol    = zero 
      ido    = 0
      info   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of PSNAUPD is used    |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | PSNAUPD.                                          |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 1
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
c        | Repeatedly call the routine PSNAUPD and take|
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call psnaupd(comm, ido, bmat, nloc, which, nev, tol, resid, 
     &        ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info )
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- OP*x                |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector, and return the matrix vector      |
c           | product to workd(ipntr(2)).               | 
c           %-------------------------------------------%
c
            call av ( comm, nloc, nx, mv_buf, 
     &                workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call PSNAUPD again.|
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
c        | Error message, check the |
c        | documentation in PSNAUPD.|
c        %--------------------------%
c
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd'
            print *, ' '
         endif
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using PSNEUPD.               |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    |
c        %-------------------------------------------%
c
         rvec = .true.
c
         call psneupd ( comm, rvec, 'A', select, d, d(1,2), v, ldv, 
     &        sigmar, sigmai, workev, bmat, nloc, which, nev, tol, 
     &        resid, ncv, v, ldv, iparam, ipntr, workd, workl,
     &        lworkl, ierr )
c
c        %-----------------------------------------------%
c        | The real part of the eigenvalue is returned   |
c        | in the first column of the two dimensional    |
c        | array D, and the imaginary part is returned   |
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
c           | Check the documentation of PSNEUPD.|
c           %------------------------------------%
c
         	if ( myid .eq. 0 ) then
             	print *, ' '
             	print *, ' Error with _neupd, info = ', ierr
             	print *, ' Check the documentation of _neupd. '
             	print *, ' '
            endif
c
         else 
c
             first  = .true.
             nconv  = iparam(5)
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
                if (d(j,2) .eq. zero)  then
c
c                  %--------------------%
c                  | Ritz value is real |
c                  %--------------------%
c
                   call av(comm, nloc, nx, mv_buf, v(1,j), ax)
                   call saxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                   d(j,3) = psnorm2( comm, nloc, ax, 1)                   
c
                else if (first) then
c
c                  %------------------------%
c                  | Ritz value is complex. |
c                  | Residual of one Ritz   |
c                  | value of the conjugate |
c                  | pair is computed.      | 
c                  %------------------------%
c        
                   call av(comm, nloc, nx, mv_buf, v(1,j), ax)
                   call saxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                   call saxpy(nloc, d(j,2), v(1,j+1), 1, ax, 1)
                   d(j,3) = psnorm2( comm, nloc, ax, 1)
                   call av(comm, nloc, nx, mv_buf, v(1,j+1), ax)
                   call saxpy(nloc, -d(j,2), v(1,j), 1, ax, 1)
                   call saxpy(nloc, -d(j,1), v(1,j+1), 1, ax, 1)
                   d(j,3) = slapy2(d(j,3), psnorm2(comm,nloc,ax,1) )
                   d(j+1,3) = d(j,3)
                   first = .false.
                else
                   first = .true.
                end if
c
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call psmout(comm, 6, nconv, 3, d, maxncv, -6,
     &            'Ritz values (Real,Imag) and direct residuals')
          end if
c
c        %-------------------------------------------%
c        | Print additional convergence information. |
c        %-------------------------------------------%
c
         if (myid .eq. 0)then
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit
     &                  Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
c
         print *, ' '
         print *, '_NDRV1 '
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
c     %---------------------------%
c     | Done with program pdndrv1.|
c     %---------------------------%
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
c     parallel matrix vector subroutine
c
c     The matrix used is the 2 dimensional convection-diffusion 
c     operator discretized using central difference.
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
c     derived from the standard central difference discretization 
c     of the 2 dimensional convection-diffusion operator 
c     (Laplacian u) + rho*(du/dx) on a unit square with zero boundary 
c     condition.
c
c     When rho*h/2 <= 1, the discrete convection-diffusion operator 
c     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX 
c     eigenvalues.
c
c     The subroutine TV is called to compute y<---T*x.
c
c----------------------------------------------------------------------------
      subroutine av (comm, nloc, nx, mv_buf, v, w)
c
c     .. MPI Declarations ...
      include           'mpif.h'
      integer           comm, nprocs, myid, ierr,
     &                  status(MPI_STATUS_SIZE)
c
      integer           nloc, nx, np, j, lo, next, prev
      Real
     &                  v(nloc), w(nloc), mv_buf(nx), one
      parameter         (one = 1.0 )
      external          saxpy, tv
c
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
c
      np = nloc/nx
      call tv(nx,v(1),w(1))
      call saxpy(nx, -one, v(nx+1), 1, w(1), 1)
c
      do 10 j = 2, np-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call saxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
         call saxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
c
      lo = (np-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call saxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
c
      next = myid + 1
      prev = myid - 1 
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v((np-1)*nx+1), nx, MPI_REAL,
     &                  next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, nx, MPI_REAL, prev, myid,
     &                  comm, status, ierr )
         call saxpy( nx, -one, mv_buf, 1, w(1), 1 )		
      endif
c
      if ( myid .gt. 0 ) then
         call mpi_send( v(1), nx, MPI_REAL,
     &                  prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, nx, MPI_REAL, next, myid,
     &                  comm, status, ierr )
         call saxpy( nx, -one, mv_buf, 1, w(lo+1), 1 )		
      endif
c
      return
      end
c=========================================================================
      subroutine tv (nx, x, y)
c
      integer           nx, j 
      Real
     &                  x(nx), y(nx), h, dd, dl, du
c
      Real
     &                  one, zero, rho
      parameter         (one = 1.0, zero = 0.0, 
     &                   rho = 0.0)
c
c     Compute the matrix vector multiplication y<---T*x
c     where T is a nx by nx tridiagonal matrix with DD on the 
c     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
c
c     When rho*h/2 <= 1, the discrete convection-diffusion operator 
c     has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX 
c     eigenvalues.
c
      h   = one / dble(nx+1)
      dd  = 4.0*one
      dl  = -one - 0.5*rho*h
      du  = -one + 0.5*rho*h
c 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end
