      program pzndrv1
!
!     Message Passing Layer: MPI
!
!     Example program to illustrate the idea of reverse communication
!     for a standard complex nonsymmetric eigenvalue problem.
!
!     We implement example one of ex-complex.doc in DOCUMENTS directory
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is obtained from the standard central difference
!         discretization of the convection-diffusion operator
!                 (Laplacian u) + rho*(du / dx)
!         in the domain omega = (0,1)x(0,1), with
!         u = 0 on the boundary of omega.
!
!     ... OP = A  and  B = I.
!     ... Assume "call av (comm, nloc, nx, mv_buf, x, y)" computes y = A*x
!     ... Use mode 1 of PZNAUPD.
!
!\BeginLib
!
!\Routines called
!     pznaupd  Parallel ARPACK reverse communication interface routine.
!     pzneupd  Parallel ARPACK routine that returns Ritz values and (optionally)
!              Ritz vectors.
!     pdznorm2 Parallel version of Level 1 BLAS that computes the norm of a complex vector.
!     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     av       Distributed matrix vector multiplication routine that computes A*x.
!     tv       Matrix vector multiplication routine that computes T*x,
!              where T is a tridiagonal matrix.  It is used in routine
!              av.
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
!     Starting Point: Complex Code FILE: ndrv1.F   SID: 2.1
!
!\SCCS Information:
! FILE: ndrv1.F   SID: 1.1   DATE OF SID: 8/13/96   RELEASE: 1
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
!
#ifdef HAVE_MPI_ICB
      use :: mpi_f08
#else
#include "mpif.h"
#endif
#include "debugF90.h"
#include "statF90.h"
!
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
      integer*4         myid, nprocs, rc, ierr, nx

!     %-----------------------------%
!     | Define maximum dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Complex*16        ax(maxn), d(maxncv),&
                        v(ldv,maxncv), workd(3*maxn),&
                        workev(3*maxncv), resid(maxn),&
                        workl(3*maxncv*maxncv+5*maxncv)
      Double precision rwork(maxncv), rd(maxncv,3)
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
                        nconv, maxitr, ishfts, mode
      integer*4         nloc
      Complex*16        sigma
      Double precision  tol
      logical           rvec
!
!     %----------------------------------------------%
!     | Local Buffers needed for MPI communication |
!     %----------------------------------------------%
!
      Complex*16        mv_buf(maxn)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision  pdznorm2
      external          pdznorm2, zaxpy
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
      mcaupd = 1
!
!     %--------------------------------------------------%
!     | The number NX is the number of interior points   |
!     | in the discretization of the 2-dimensional       |
!     | convection-diffusion operator on the unit        |
!     | square with zero Dirichlet boundary condition.   |
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%
!
      nx    = 10
      n     = nx*nx
      nev   = 4
      ncv   = 20
!
!     %--------------------------------------%
!     | Set up distribution of data to nodes |
!     %--------------------------------------%
!
      nloc = (nx / nprocs)*nx
      if ( mod(nx, nprocs) .gt. myid ) nloc = nloc + nx
!
      if ( nloc .gt. maxn ) then
         print *, ' ERROR with _NDRV1: NLOC is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'LM'
!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD as         |
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         |
!     %---------------------------------------------------%
!
      lworkl  = 3*ncv**2+5*ncv
      tol    = 0.0
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 1
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
!        | Repeatedly call the routine ZNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call pznaupd ( comm, ido, bmat, nloc, which,&
              nev, tol, resid, ncv, v, ldv, iparam, ipntr,&
              workd, workl, lworkl, rwork,info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               |
!           %-------------------------------------------%
!
            call av ( comm, nloc, nx, mv_buf,&
                      workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
         end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD  |
!        %--------------------------%
!
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd'
            print *, ' '
         endif
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call pzneupd (comm, rvec, 'A', select, d, v, ldv, sigma,&
              workev, bmat, nloc, which, nev, tol, resid, ncv,&
              v, ldv, iparam, ipntr, workd, workl, lworkl,&
              rwork, info)
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( info .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD. |
!           %------------------------------------%
!
            if ( myid .eq. 0 ) then
                print *, ' '
                print *, ' Error with _neupd, info = ', info
                print *, ' Check the documentation of _neupd. '
                print *, ' '
            endif
!
         else
!
             nconv = iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                call av(comm, nloc, nx, mv_buf, v(1,j), ax)
                call zaxpy(nloc, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble(d(j))
                rd(j,2) = dimag(d(j))
                rd(j,3) = pdznorm2(comm, nloc, ax, 1)
!
 20          continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
             call pdmout(comm, 6, nconv, 3, rd, maxncv, -6,&
                  'Ritz values (Real, Imag) and direct residuals')
          end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
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
         print *, '_NDRV1'
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
!     | Done with program pzndrv1. |
!     %----------------------------%
!
 9000 continue
!
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
!     parallel matrix vector subroutine
!
!     The matrix used is the convection-diffusion operator
!     discretized using centered difference.
!
!     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
!     tridiagonal matrix
!
!                  | T -I          |
!                  |-I  T -I       |
!             OP = |   -I  T       |
!                  |        ...  -I|
!                  |           -I T|
!
!     derived from the standard central difference  discretization
!     of the convection-diffusion operator (Laplacian u) + rho*(du/dx)
!     with zero boundary condition.
!
!     The subroutine TV is called to computed y<---T*x.
!
!----------------------------------------------------------------------------
      subroutine av (comm, nloc, nx, mv_buf, v, w )
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
      integer*4         nprocs, myid, ierr, next, prev, tag
!
      integer*4         nloc, nx, np
      integer           j, lo
      Complex*16        v(nloc), w(nloc), mv_buf(nx), one
      parameter         (one = (1.0, 0.0))
      external          zaxpy, tv
!
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )
!
      np = nloc/nx
      call tv(nx,v(1),w(1))
      call zaxpy(nx, -one, v(nx+1), 1, w(1), 1)
!
      do 10 j = 2, np-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call zaxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
         call zaxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue
!
      lo = (np-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call zaxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
!
      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         tag = myid+1
         call mpi_send( v((np-1)*nx+1), nx, MPI_DOUBLE_COMPLEX,&
                        next, tag, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, nx, MPI_DOUBLE_COMPLEX, prev, myid,&
                        comm, status, ierr )
         call zaxpy( nx, -one, mv_buf, 1, w(1), 1 )
      endif
!
      if ( myid .gt. 0 ) then
         tag = myid-1
         call mpi_send( v(1), nx, MPI_DOUBLE_COMPLEX,&
                        prev, tag, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, nx, MPI_DOUBLE_COMPLEX, next, myid,&
                        comm, status, ierr )
         call zaxpy( nx, -one, mv_buf, 1, w(lo+1), 1 )
      endif
!
      return
      end
!=========================================================================
      subroutine tv (nx, x, y)
!
      integer*4         nx
      integer           j
      Complex*16        x(nx), y(nx), h, dd, dl, du
!
      Complex*16        one, rho
      parameter         (one = (1.0, 0.0), rho = (100.0, 0.0))
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal
!
      h   = one / dcmplx(nx+1)
      dd  = (4.0, 0.0)
      dl  = -one - (0.5, 0.0)*rho*h
      du  = -one + (0.5, 0.0)*rho*h
!
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
 10   continue
      y(nx) =  dl*x(nx-1) + dd*x(nx)
      return
      end
