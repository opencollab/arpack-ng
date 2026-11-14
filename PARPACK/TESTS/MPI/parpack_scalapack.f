      program parpack_scalapack
c
c     PARPACK and ScaLAPACK linked into the same program.
c
c     ScaLAPACK provides routines whose names were also used by
c     PARPACK/ARPACK with a different interface:
c       - PDLAMCH/PSLAMCH (BLACS context argument), PARPACK versions
c         take an MPI communicator: renamed PDLAMCH10/PSLAMCH10 (#21),
c       - ZZDOTC/CCDOTC (subroutines used by PZDOTC/PCDOTC), ARPACK
c         versions are functions: renamed ZZDOTC10/CCDOTC10 (#505).
c     With a name collision, depending on the link order, either the
c     PARPACK solver or the ScaLAPACK routine calls the wrong code.
c     To make such a mix-up visible:
c       - this test is linked with ScaLAPACK before PARPACK,
c       - PARPACK uses its own communicator (MPI_COMM_DUP), so its
c         handle differs from the BLACS context.
c
c     1. Check ScaLAPACK p[sd]lamch against LAPACK [sd]lamch.
c     2. For both complex precisions (p[cz] routines use p[sd]lamch10
c        and [cz][cz]dotc10):
c        a. compute the 2 largest eigenpairs of A = diag(1,...,24) with
c           p[cz]naupd/p[cz]neupd (expected eigenvalues: 24 and 23),
c        b. check the computed eigenvectors with PBLAS p[cz]dotc using
c           a ScaLAPACK descriptor matching the PARPACK distribution:
c           <z1,z1> = 1, <z1,z2> = 0 and <z1,A*z1> = lambda1.
c
c     Run with 1, 2, 3 or 4 MPI processes (NCV must not exceed local N).
c
      implicit none
      include 'mpif.h'
      integer          nglob, nev, ncv, ldv, lworkl
      parameter        (nglob = 24, nev = 2, ncv = 5, ldv = nglob,
     &                  lworkl = 3*ncv*ncv+5*ncv)
      integer          iparam(11), ipntr(14), desc(9), ido, info, ierr,
     &                 comm, myid, nprocs, n, i, pass, row, ictxt,
     &                 nprow, npcol, myrow, mycol
      logical          select(ncv), failed
      complex          resid(ldv), v(ldv,ncv), workd(3*ldv),
     &                 workl(lworkl), workev(2*ncv), d(nev+1),
     &                 w(ldv), cdot
      real             rwork(ncv), ctol
      double complex   zresid(ldv), zv(ldv,ncv), zworkd(3*ldv),
     &                 zworkl(lworkl), zworkev(2*ncv), zd(nev+1),
     &                 zw(ldv), zdot
      double precision zrwork(ncv), dlam(2), tol, ztol
      real             pslamch, slamch
      double precision pdlamch, dlamch
      external         pslamch, slamch, pdlamch, dlamch
c
      call MPI_INIT(ierr)
      call MPI_COMM_DUP(MPI_COMM_WORLD, comm, ierr)
      call MPI_COMM_RANK(comm, myid, ierr)
      call MPI_COMM_SIZE(comm, nprocs, ierr)
      if (mod(nglob,nprocs) .ne. 0) call MPI_ABORT(comm, 2, ierr)
      n = nglob / nprocs
c
c     ScaLAPACK view of the PARPACK distribution: nprocs x 1 process
c     grid, process p owns rows p*n+1 ... p*n+n of each vector.
c
      call BLACS_GET(-1, 0, ictxt)
      call BLACS_GRIDINIT(ictxt, 'Row-major', nprocs, 1)
      call BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
      if (myrow .ne. myid) call MPI_ABORT(comm, 2, ierr)
      call DESCINIT(desc, nglob, 1, n, 1, 0, 0, ictxt, ldv, info)
      if (info .ne. 0) call MPI_ABORT(comm, 2, ierr)
c
      failed = .false.
c
c     ScaLAPACK machine parameters (BLACS context).
c
      if (pslamch(ictxt, 'E') .ne. slamch('E') .or.
     &    pdlamch(ictxt, 'E') .ne. dlamch('E')) then
         print *, ' Error: ScaLAPACK p[sd]lamch do not match [sd]lamch'
         failed = .true.
      end if
c
      do 100 pass = 1, 2
         ido = 0
         info = 1
c        tol <= 0: PARPACK uses machine precision (and updates tol).
         ctol = 0.0E+0
         ztol = 0.0D+0
         do 10 i = 1, n
            resid(i) = cmplx(1.0E+0, 0.0E+0)
            zresid(i) = dcmplx(1.0D+0, 0.0D+0)
   10    continue
         do 20 i = 1, 11
            iparam(i) = 0
   20    continue
         iparam(1) = 1
         iparam(3) = 300
         iparam(7) = 1
c
   30    continue
         if (pass .eq. 1) then
            call pcnaupd(comm, ido, 'I', n, 'LM', nev, ctol, resid,
     &                   ncv, v, ldv, iparam, ipntr, workd, workl,
     &                   lworkl, rwork, info)
         else
            call pznaupd(comm, ido, 'I', n, 'LM', nev, ztol, zresid,
     &                   ncv, zv, ldv, iparam, ipntr, zworkd, zworkl,
     &                   lworkl, zrwork, info)
         end if
         if (ido .eq. -1 .or. ido .eq. 1) then
            do 40 i = 1, n
               row = myid*n+i
               if (pass .eq. 1) then
                  workd(ipntr(2)+i-1) = real(row)*workd(ipntr(1)+i-1)
               else
                  zworkd(ipntr(2)+i-1) = dble(row)*zworkd(ipntr(1)+i-1)
               end if
   40       continue
            go to 30
         end if
         if (ido .ne. 99 .or. info .ne. 0 .or.
     &       iparam(5) .ne. nev) then
            print *, ' Error with p[cz]naupd, info = ', info,
     &               ', nconv = ', iparam(5)
            call MPI_ABORT(comm, 2, ierr)
         end if
c
         if (pass .eq. 1) then
            call pcneupd(comm, .true., 'A', select, d, v, ldv,
     &                   (0.0E+0,0.0E+0), workev, 'I', n, 'LM', nev,
     &                   ctol, resid, ncv, v, ldv, iparam, ipntr,
     &                   workd, workl, lworkl, rwork, info)
         else
            call pzneupd(comm, .true., 'A', select, zd, zv, ldv,
     &                   (0.0D+0,0.0D+0), zworkev, 'I', n, 'LM', nev,
     &                   ztol, zresid, ncv, zv, ldv, iparam, ipntr,
     &                   zworkd, zworkl, lworkl, zrwork, info)
         end if
         if (info .ne. 0) then
            print *, ' Error with p[cz]neupd, info = ', info
            call MPI_ABORT(comm, 2, ierr)
         end if
c
c        Eigenvalues: {24, 23} in any order, imaginary parts zero.
c
         if (pass .eq. 1) then
            tol = 1.0D-3
            dlam(1) = dble(real(d(1)))
            dlam(2) = dble(real(d(2)))
            if (abs(aimag(d(1))) .gt. tol .or.
     &          abs(aimag(d(2))) .gt. tol) failed = .true.
         else
            tol = 1.0D-10
            dlam(1) = dble(zd(1))
            dlam(2) = dble(zd(2))
            if (abs(dimag(zd(1))) .gt. tol .or.
     &          abs(dimag(zd(2))) .gt. tol) failed = .true.
         end if
         if (abs(max(dlam(1),dlam(2)) - dble(nglob)) .gt. tol .or.
     &       abs(min(dlam(1),dlam(2)) - dble(nglob-1)) .gt. tol)
     &      failed = .true.
         if (myid .eq. 0) print *, ' Pass ', pass, ' eigenvalues: ',
     &                             dlam(1), dlam(2)
c
c        Check the eigenvectors with ScaLAPACK (PBLAS).
c
         do 50 i = 1, n
            row = myid*n+i
            if (pass .eq. 1) then
               w(i) = real(row)*v(i,1)
            else
               zw(i) = dble(row)*zv(i,1)
            end if
   50    continue
c
         if (pass .eq. 1) then
            call pcdotc(nglob, cdot, v(1,1), 1, 1, desc, 1,
     &                  v(1,1), 1, 1, desc, 1)
            zdot = dcmplx(dble(real(cdot)), dble(aimag(cdot)))
         else
            call pzdotc(nglob, zdot, zv(1,1), 1, 1, desc, 1,
     &                  zv(1,1), 1, 1, desc, 1)
         end if
         if (myid .eq. 0) print *, ' Pass ', pass, ' <z1,z1>   = ', zdot
         if (abs(zdot - dcmplx(1.0D+0,0.0D+0)) .gt. tol)
     &      failed = .true.
c
         if (pass .eq. 1) then
            call pcdotc(nglob, cdot, v(1,1), 1, 1, desc, 1,
     &                  v(1,2), 1, 1, desc, 1)
            zdot = dcmplx(dble(real(cdot)), dble(aimag(cdot)))
         else
            call pzdotc(nglob, zdot, zv(1,1), 1, 1, desc, 1,
     &                  zv(1,2), 1, 1, desc, 1)
         end if
         if (myid .eq. 0) print *, ' Pass ', pass, ' <z1,z2>   = ', zdot
         if (abs(zdot) .gt. tol) failed = .true.
c
         if (pass .eq. 1) then
            call pcdotc(nglob, cdot, v(1,1), 1, 1, desc, 1,
     &                  w, 1, 1, desc, 1)
            zdot = dcmplx(dble(real(cdot)), dble(aimag(cdot)))
         else
            call pzdotc(nglob, zdot, zv(1,1), 1, 1, desc, 1,
     &                  zw, 1, 1, desc, 1)
         end if
         if (myid .eq. 0) print *, ' Pass ', pass, ' <z1,Az1>  = ', zdot
         if (abs(zdot - dcmplx(dlam(1),0.0D+0)) .gt. tol*dble(nglob))
     &      failed = .true.
  100 continue
c
      call BLACS_GRIDEXIT(ictxt)
      call BLACS_EXIT(1)
      call MPI_COMM_FREE(comm, ierr)
      call MPI_FINALIZE(ierr)
      if (failed) then
         print *, ' Error: PARPACK and ScaLAPACK results do not match'
         stop 1
      end if
      end
