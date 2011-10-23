      program psntest1 
c
c     Message Passing Layer: BLACS
c
c     Example program to illustrate the idea of reverse communication
c     for a standard nonsymmetric eigenvalue problem.
c
c     We implement example one of ex-nonsym.doc in DOCUMENTS directory
c
c\Test-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is random diagonal matrix with 4 separated eigenvalues.
c     ... OP = A  and  B = I.
c     ... Assume "call av ( nloc, diag, x, y)" computes y = A*x.
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
c
c\Author
c     Kristi Maschhoff 
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: 
c FILE: %M%   SID: %I%   DATE OF SID: %G%   RELEASE: %R%
c
c\Remarks
c     1. None
c
c\EndLib
c---------------------------------------------------------------------------
c
      include 'debug.h'
      include 'stat.h'
 
c     %-----------------%
c     | BLACS INTERFACE |
c     %-----------------%
c
      integer           comm, iam, nprocs, nloc,
     &                  nprow, npcol, myprow, mypcol
c
      external          BLACS_PINFO, BLACS_SETUP, BLACS_GET,
     &                  BLACS_GRIDINIT, BLACS_GRIDINFO
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
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14), iseed(4)
      logical           select(maxncv)
      Real
     &                  ax(maxn), d(maxncv,3), resid(maxn), diag(maxn), 
     &                  v(ldv,maxncv), workd(3*maxn), 
     &                  workev(3*maxncv), 
     &                  workl(3*maxncv*maxncv+6*maxncv)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode, idist
      Real
     &                  tol, sigmar, sigmai
      logical           first, rvec
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &                  zero, one
      parameter         (zero = 0.0, one = 1.0)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real
     &                  slapy2, psnorm2
      external          slapy2, saxpy, psnorm2, slarnv
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
      call BLACS_PINFO( iam, nprocs )
c
c     If in PVM, create virtual machine if it doesn't exist
c
      if (nprocs .lt. 1) then
         if (iam .eq. 0) then
              write(*,1000)
              read(*, 2000) nprocs
         endif
         call BLACS_SETUP( iam, nprocs )
      endif
c
1000  format('How many processes in machine?')
2000  format(I3)
c
c     Set up processors in 1D Grid
c
      nprow = nprocs
      npcol = 1
c
c     Get default system context, and define grid
c
      call BLACS_GET( 0, 0, comm )
      call BLACS_GRIDINIT( comm, 'Row', nprow, npcol )
      call BLACS_GRIDINFO( comm, nprow, npcol, myprow, mypcol )
c
c     If I'm not in grid, go to end of program
c
      if ( (myprow .ge. nprow) .or. (mypcol .ge. npcol) ) goto 9000
c
      ndigit = -3
      logfil = 6
      mnaupd = 1
c
      n     = 10000*nprocs
      nev   = 4
      ncv   = 20 
c
c     %--------------------------------------%
c     | Set up distribution of data to nodes |
c     %--------------------------------------%
c
      nloc = 10000
c
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
      which = 'SM'
c
c     %-----------------------------------%
c     | Generate random diagonal matrix   |
c     | Isolate 4 extreamal eigenvalues   |
c     %-----------------------------------%
c
      idist = 1
      iseed(1) = 15
      iseed(2) = 35
      iseed(3) = 52
      iseed(4) = 7
      call slarnv ( idist, iseed, nloc, diag )
      diag(1) = diag(1) + 1.01
      diag(2) = diag(2) + 1.01
      diag(3) = diag(3) + 1.01
      diag(4) = diag(4) + 1.01
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
      info   = 1
      do 50 j=1,nloc
         resid(j) = one
 50   continue
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
            call av ( nloc, diag, workd(ipntr(1)), workd(ipntr(2)))
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
         if ( myprow .eq. 0 ) then
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
         	if ( myprow .eq. 0 ) then
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
                   call av( nloc, diag, v(1,j), ax)
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
                   call av( nloc, diag, v(1,j), ax)
                   call saxpy(nloc, -d(j,1), v(1,j), 1, ax, 1)
                   call saxpy(nloc, d(j,2), v(1,j+1), 1, ax, 1)
                   d(j,3) = psnorm2( comm, nloc, ax, 1)
                   call av( nloc, diag, v(1,j+1), ax)
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
         if (myprow .eq. 0)then
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
c     | Release resources BLACS |
c     %-------------------------%
c
      call BLACS_GRIDEXIT ( comm )
      call BLACS_EXIT(0)
c
      end
c 
c==========================================================================
c
c     parallel matrix vector subroutine
c
      subroutine av (n, diag, v, w)
      integer           n, j
      Real
     &                  v(n), w(n), diag(n)
c
      do 10 j = 1, n
         w(j) = diag(j)*v(j)
  10  continue
c
      return
      end
