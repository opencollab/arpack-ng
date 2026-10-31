      program parpack_default_tol
c
c     Check the plain Fortran default TOL across reverse communication
c     on the same matrices as icb_parpack_default_tol.c.
c
      implicit none
      include 'mpif.h'
      integer          comm, myid, nprocs, n, start, ierr
      call MPI_INIT(ierr)
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm, myid, ierr)
      call MPI_COMM_SIZE(comm, nprocs, ierr)
      if (nprocs .gt. 10) call MPI_ABORT(comm, 1, ierr)
      n = 10 / nprocs
      if (myid .lt. mod(10,nprocs)) n = n + 1
      start = myid*(10/nprocs) + min(myid,mod(10,nprocs))
      call solve(comm, myid, nprocs, n, start, .false., 0.0D+0)
      call solve(comm, myid, nprocs, n, start, .false., -1.0D+0)
      call solve(comm, myid, nprocs, n, start, .true., 0.0D+0)
      call solve(comm, myid, nprocs, n, start, .true., -1.0D+0)
      call MPI_FINALIZE(ierr)
      end
c
      subroutine solve(comm, myid, nprocs, n, start, upper, initial_tol)
      implicit none
      include 'mpif.h'
      logical          upper, select(6)
      integer          ntotal, nev, ncv, lworkl
      parameter        (ntotal = 10, nev = 2, ncv = 6,
     &                  lworkl = 3*ncv*(ncv+2))
      integer          comm, myid, nprocs, n, start, ierr,
     &                 ido, info, iparam(11), ipntr(14), i, j
      double precision initial_tol, tol, eps, pdlamch10, expected,
     &                 resid(ntotal), v(ntotal,ncv), workd(3*ntotal),
     &                 workl(lworkl), workev(3*ncv), dr(nev+1),
     &                 di(nev+1), z(ntotal,nev+1), y(ntotal),
     &                 sums(2), totals(2), r
      external         pdlamch10
c
c     A has diagonal 1,...,10 and, if UPPER, superdiagonal 1.
c     Use a deterministic starting vector distributed over the ranks.
c
      ido = 0
      info = 1
      tol = initial_tol
      eps = pdlamch10(comm, 'EpsMach')
      do 10 i = 1, n
         resid(i) = dble(start+i)
   10 continue
      do 20 i = 1, 11
         iparam(i) = 0
   20 continue
      iparam(1) = 1
      iparam(3) = 100
      iparam(4) = 1
      iparam(7) = 1
c
   30 continue
      if (upper) then
         call pdnaupd(comm, ido, 'I', n, 'LM', nev, tol, resid, ncv,
     &                v, ntotal, iparam, ipntr, workd, workl,
     &                lworkl, info)
      else
         call pdsaupd(comm, ido, 'I', n, 'LM', nev, tol, resid, ncv,
     &                v, ntotal, iparam, ipntr, workd, workl,
     &                lworkl, info)
      end if
      if (tol .ne. eps) then
         print *, ' AUPD did not retain default TOL: ', tol
         call MPI_ABORT(comm, 1, ierr)
      end if
      if (ido .eq. -1 .or. ido .eq. 1) then
         call matvec(comm, myid, nprocs, n, start, upper,
     &               workd(ipntr(1)), workd(ipntr(2)))
         go to 30
      end if
      if (ido .ne. 99 .or. info .ne. 0 .or.
     &    iparam(5) .ne. nev) then
         print *, ' AUPD: upper, initial TOL, IDO, INFO, NCONV = ',
     &            upper, initial_tol, ido, info, iparam(5)
         call MPI_ABORT(comm, 1, ierr)
      end if
c
c     Pass the same TOL variable to EUPD, as in the Fortran examples.
c
      if (upper) then
         call pdneupd(comm, .true., 'A', select, dr, di, z, ntotal,
     &                0.0D+0, 0.0D+0, workev, 'I', n, 'LM', nev, tol,
     &                resid, ncv, v, ntotal, iparam, ipntr, workd,
     &                workl, lworkl, info)
      else
         call pdseupd(comm, .true., 'A', select, dr, z, ntotal,
     &                0.0D+0, 'I', n, 'LM', nev, tol, resid, ncv, v,
     &                ntotal, iparam, ipntr, workd, workl, lworkl, info)
         di(1) = 0.0D+0
         di(2) = 0.0D+0
      end if
      if (info .ne. 0 .or. tol .ne. eps) then
         print *, ' EUPD: INFO, TOL = ', info, tol
         call MPI_ABORT(comm, 1, ierr)
      end if
c
c     Check eigenvalues and globally reduced norms and residuals with
c     the same bounds as the ICB test, including the boundary rows.
c
      do 50 j = 1, nev
         expected = dble(ntotal-nev+j)
         if (upper .and. dr(1) .gt. dr(2)) then
            expected = dble(ntotal+1-j)
         end if
         call matvec(comm, myid, nprocs, n, start, upper, z(1,j), y)
         sums(1) = 0.0D+0
         sums(2) = 0.0D+0
         do 40 i = 1, n
            r = y(i) - dr(j)*z(i,j)
            sums(1) = sums(1) + r*r
            sums(2) = sums(2) + z(i,j)*z(i,j)
   40    continue
         call MPI_ALLREDUCE(sums, totals, 2, MPI_DOUBLE_PRECISION,
     &                      MPI_SUM, comm, ierr)
         if (.not. (abs(dr(j)-expected) .le. 1.0D-10 .and.
     &              abs(di(j)) .le. 1.0D-10 .and.
     &              totals(2) .ge. 0.9D+0 .and.
     &              totals(1) .le. 1.0D-20*totals(2))) then
            print *, ' Incorrect eigenpair: ', j, dr(j), di(j), totals
            call MPI_ABORT(comm, 1, ierr)
         end if
   50 continue
      if (myid .eq. 0) then
         print *, ' upper, initial TOL, Ritz values: ', upper,
     &            initial_tol, dr(1), dr(2)
      end if
      end
c
      subroutine matvec(comm, myid, nprocs, n, start, upper, x, y)
      implicit none
      include 'mpif.h'
      logical          upper
      integer          comm, myid, nprocs, n, start, i, ierr,
     &                 prev, next, status(MPI_STATUS_SIZE)
      double precision x(n), y(n), halo
c
c     The last local row needs the first entry from the next rank.
c
      halo = 0.0D+0
      if (upper) then
         prev = MPI_PROC_NULL
         next = MPI_PROC_NULL
         if (myid .gt. 0) prev = myid - 1
         if (myid+1 .lt. nprocs) next = myid + 1
         call MPI_SENDRECV(x, 1, MPI_DOUBLE_PRECISION, prev, 0,
     &                     halo, 1, MPI_DOUBLE_PRECISION, next, 0,
     &                     comm, status, ierr)
      end if
      do 10 i = 1, n
         y(i) = dble(start+i)*x(i)
         if (upper) then
            if (i .lt. n) then
               y(i) = y(i) + x(i+1)
            else
               y(i) = y(i) + halo
            end if
         end if
   10 continue
      end
