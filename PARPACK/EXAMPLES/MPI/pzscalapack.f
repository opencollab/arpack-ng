      program pzscalapack
c
c     Example: PARPACK together with ScaLAPACK.
c
c     The operator is a dense complex N x N matrix A stored and
c     distributed with ScaLAPACK, and the matrix-vector product
c     y <- A*x needed by PARPACK is computed with PBLAS PZGEMV.
c     The computed eigenpairs are then checked with PBLAS:
c     ||A*z - lambda*z|| / |lambda| with PZGEMV, PZAXPY and PDZNRM2.
c
c     A is the non-Hermitian tridiagonal Toeplitz matrix
c        A(i,i) = 2, A(i,i+1) = -1 + i/2, A(i+1,i) = -1
c     whose eigenvalues are
c        2 + 2*sqrt((-1+i/2)*(-1)) * cos(k*pi/(N+1)), k = 1, ..., N.
c
c     Data distribution
c     -----------------
c     PARPACK distributes vectors by contiguous blocks of rows, one
c     block per MPI process. The same layout is described for ScaLAPACK
c     by a NPROCS x 1 process grid and a row block size MB:
c        - process p owns global rows p*MB+1 ... min((p+1)*MB, N),
c        - A uses descriptor (N, N, MB, N): each process stores whole
c          rows of A,
c        - a vector uses descriptor (N, 1, MB, 1), so the local part
c          of a PARPACK vector can be passed to PBLAS as is.
c
c     Routines called
c     ---------------
c     pznaupd  PARPACK reverse communication interface routine.
c     pzneupd  PARPACK routine that returns Ritz values and vectors.
c     pzgemv   PBLAS distributed matrix-vector product.
c     pzaxpy   PBLAS distributed y <- alpha*x + y.
c     pdznrm2  PBLAS distributed 2-norm.
c     descinit ScaLAPACK descriptor initialization.
c     numroc   ScaLAPACK number of local rows.
c     blacs_*  BLACS process grid routines.
c
c     Run with any number of MPI processes such that NCV does not exceed
c     the number of local rows, e.g. "mpirun -n 2 pzscalapack".
c
      implicit none
      include 'mpif.h'
c
      integer          maxn, maxnev, maxncv, ldv, lworkl
      parameter        (maxn = 256, maxnev = 10, maxncv = 30,
     &                  ldv = maxn, lworkl = 3*maxncv*maxncv+5*maxncv)
c
      integer          iparam(11), ipntr(14), desca(9), descx(9),
     &                 comm, myid, nprocs, ierr, ictxt, nprow, npcol,
     &                 myrow, mycol, n, nloc, mb, nev, ncv, ido, info,
     &                 i, j, gi, nconv, maxitr, mode, numroc
      logical          select(maxncv), rvec, failed
      character        bmat*1, which*2
      Double complex   a(ldv,maxn), v(ldv,maxncv), workd(3*ldv),
     &                 workl(lworkl), workev(2*maxncv), d(maxnev+1),
     &                 resid(ldv), ax(ldv), sigma, diag, sup, sub,
     &                 one, zero
      Double precision rwork(maxncv), tol, rnorm
      parameter        (one = (1.0D+0,0.0D+0),
     &                  zero = (0.0D+0,0.0D+0))
      external         numroc
c
      call MPI_INIT(ierr)
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm, myid, ierr)
      call MPI_COMM_SIZE(comm, nprocs, ierr)
c
c     Problem size and PARPACK parameters.
c
      n     = 100
      nev   = 4
      ncv   = 20
      bmat  = 'I'
      which = 'LM'
c
c     NPROCS x 1 BLACS grid: grid row p is MPI rank p.
c
      call BLACS_GET(-1, 0, ictxt)
      call BLACS_GRIDINIT(ictxt, 'Row-major', nprocs, 1)
      call BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
c
      mb   = (n + nprocs - 1) / nprocs
      nloc = numroc(n, mb, myrow, 0, nprow)
      if (nloc .lt. ncv .or. nloc .gt. ldv) then
         if (myid .eq. 0) print *, ' ERROR: need NCV <= local N <= LDV'
         call MPI_ABORT(comm, 1, ierr)
      end if
c
      call descinit(desca, n, n, mb, n, 0, 0, ictxt, ldv, info)
      call descinit(descx, n, 1, mb, 1, 0, 0, ictxt, ldv, info)
c
c     Fill the local rows of A.
c
      diag = ( 2.0D+0, 0.0D+0)
      sup  = (-1.0D+0, 0.5D+0)
      sub  = (-1.0D+0, 0.0D+0)
      do 20 j = 1, n
         do 10 i = 1, nloc
            a(i,j) = zero
   10    continue
   20 continue
      do 30 i = 1, nloc
         gi = myrow*mb + i
         a(i,gi) = diag
         if (gi .lt. n) a(i,gi+1) = sup
         if (gi .gt. 1) a(i,gi-1) = sub
   30 continue
c
c     PARPACK reverse communication loop.
c
      tol    = 0.0D+0
      sigma  = zero
      ido    = 0
      info   = 0
      maxitr = 300
      mode   = 1
      do 40 i = 1, 11
         iparam(i) = 0
   40 continue
      iparam(1) = 1
      iparam(3) = maxitr
      iparam(7) = mode
c
   50 continue
         call pznaupd(comm, ido, bmat, nloc, which, nev, tol, resid,
     &                ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,
     &                rwork, info)
         if (ido .eq. -1 .or. ido .eq. 1) then
c           y <- A*x with x = workd(ipntr(1)), y = workd(ipntr(2)).
            call pzgemv('N', n, n, one, a, 1, 1, desca,
     &                  workd(ipntr(1)), 1, 1, descx, 1, zero,
     &                  workd(ipntr(2)), 1, 1, descx, 1)
            go to 50
         end if
c
      if (info .lt. 0) then
         if (myid .eq. 0) print *, ' Error with pznaupd, info = ', info
         call MPI_ABORT(comm, 1, ierr)
      end if
c
      rvec = .true.
      call pzneupd(comm, rvec, 'A', select, d, v, ldv, sigma, workev,
     &             bmat, nloc, which, nev, tol, resid, ncv, v, ldv,
     &             iparam, ipntr, workd, workl, lworkl, rwork, ierr)
      if (ierr .ne. 0) then
         if (myid .eq. 0) print *, ' Error with pzneupd, info = ', ierr
         call MPI_ABORT(comm, 1, ierr)
      end if
c
c     Check each Ritz pair with PBLAS: ||A*z - lambda*z|| / |lambda|.
c
      nconv  = iparam(5)
      failed = .false.
      if (myid .eq. 0) then
         print *, ' '
         print *, ' PARPACK + ScaLAPACK example (pzscalapack)'
         print *, ' Size of the matrix is ', n
         print *, ' Number of processes is ', nprocs
         print *, ' Number of converged Ritz values is ', nconv
         print *, ' Ritz value, relative residual:'
      end if
      do 60 j = 1, nconv
         call pzgemv('N', n, n, one, a, 1, 1, desca, v(1,j), 1, 1,
     &               descx, 1, zero, ax, 1, 1, descx, 1)
         call pzaxpy(n, -d(j), v(1,j), 1, 1, descx, 1,
     &               ax, 1, 1, descx, 1)
         call pdznrm2(n, rnorm, ax, 1, 1, descx, 1)
         rnorm = rnorm / abs(d(j))
         if (myid .eq. 0) print *, '   ', d(j), rnorm
         if (rnorm .gt. 1.0D-8) failed = .true.
   60 continue
      if (nconv .lt. nev) failed = .true.
c
      call BLACS_GRIDEXIT(ictxt)
      call BLACS_EXIT(1)
      call MPI_FINALIZE(ierr)
      if (failed) then
         print *, ' Error: unexpected Ritz pairs'
         stop 1
      end if
      end
