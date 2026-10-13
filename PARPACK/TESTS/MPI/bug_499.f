      program parpack_bug_499
c
c     Test eigenvalues-only 'BE' extraction in both PARPACK precisions.
c     A = diag(1,2,3,4), NEV = 2: the expected Ritz values are 1 and 4.
c     Distribute the four rows across one, two or four MPI processes.
c
      implicit none
      include 'mpif.h'
      integer          n, nev, ncv, ldv, lworkl
      parameter        (nev = 2, ncv = 4, ldv = 4,
     &                  lworkl = ncv*(ncv+8))
      integer          iparam(11), ipntr(11), ido, info, ierr, i,
     &                 comm, myid, nprocs, pass, row
      logical          select(ncv), failed
      real             resid(ldv), v(ldv,ncv), workd(3*ldv),
     &                 workl(lworkl), d(nev), z(ldv,nev)
      double precision dresid(ldv), dv(ldv,ncv), dworkd(3*ldv),
     &                 dworkl(lworkl), dd(nev), dz(ldv,nev)
c
      call MPI_INIT(ierr)
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm, myid, ierr)
      call MPI_COMM_SIZE(comm, nprocs, ierr)
      if (mod(ldv,nprocs) .ne. 0) call MPI_ABORT(comm, 2, ierr)
      n = ldv / nprocs
      failed = .false.
      do 50 pass = 1, 2
         ido = 0
         info = 1
         do 10 i = 1, n
            resid(i) = real(myid*n+i)
            dresid(i) = dble(myid*n+i)
   10    continue
         do 20 i = 1, 11
            iparam(i) = 0
   20    continue
         iparam(1) = 1
         iparam(3) = 100
         iparam(7) = 1
c
   30    continue
         if (pass .eq. 1) then
            call pssaupd (comm, ido, 'I', n, 'BE', nev, 1.0E-5,
     &                   resid, ncv, v, ldv, iparam, ipntr, workd,
     &                   workl, lworkl, info)
         else
            call pdsaupd (comm, ido, 'I', n, 'BE', nev, 1.0D-5,
     &                   dresid, ncv, dv, ldv, iparam, ipntr, dworkd,
     &                   dworkl, lworkl, info)
         end if
         if (ido .eq. -1 .or. ido .eq. 1) then
            do 40 i = 1, n
               row = myid*n+i
               if (pass .eq. 1) then
                  workd(ipntr(2)+i-1) = real(row)*workd(ipntr(1)+i-1)
               else
                  dworkd(ipntr(2)+i-1) = dble(row)*dworkd(ipntr(1)+i-1)
               end if
   40       continue
            go to 30
         end if
         if (ido .ne. 99 .or. info .ne. 0 .or.
     &       iparam(5) .ne. nev) then
            print *, ' Error with _saupd, info = ', info
            call MPI_ABORT(comm, 2, ierr)
         end if
c
         if (pass .eq. 1) then
            call psseupd (comm, .false., 'A', select, d, z, ldv,
     &                   0.0E+0, 'I', n, 'BE', nev, 1.0E-5, resid,
     &                   ncv, v, ldv, iparam, ipntr, workd, workl,
     &                   lworkl, info)
            if (myid .eq. 0) print *, ' Single Ritz values: ', d
            if (.not. (abs(d(1)-1.0E+0) .le. 1.0E-4 .and.
     &                 abs(d(2)-4.0E+0) .le. 1.0E-4)) failed = .true.
         else
            call pdseupd (comm, .false., 'A', select, dd, dz, ldv,
     &                   0.0D+0, 'I', n, 'BE', nev, 1.0D-5, dresid,
     &                   ncv, dv, ldv, iparam, ipntr, dworkd, dworkl,
     &                   lworkl, info)
            if (myid .eq. 0) print *, ' Double Ritz values: ', dd
            if (.not. (abs(dd(1)-1.0D+0) .le. 1.0D-4 .and.
     &                 abs(dd(2)-4.0D+0) .le. 1.0D-4)) failed = .true.
         end if
         if (info .ne. 0) then
            print *, ' Error with _seupd, info = ', info
            call MPI_ABORT(comm, 2, ierr)
         end if
   50 continue
      call MPI_FINALIZE(ierr)
      if (failed) then
         print *, ' Error: expected Ritz values 1 and 4'
         stop 1
      end if
      end
