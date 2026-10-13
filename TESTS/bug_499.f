      program bug_499
c
c     Test for issue #499: ssaupd with WHICH = 'BE' must return Ritz
c     values from both ends of the spectrum. With NCV = N, ssaup2
c     reduced NP to zero before the 'BE' swap, so sseupd with
c     RVEC = .false. returned the two largest eigenvalues instead of
c     the smallest and the largest one.
c
c     A = diag(1,2,3,4), NEV = 2: the expected eigenvalues are 1 and 4.
c
      implicit none
      integer          n, nev, ncv, ldv, lworkl
      parameter        (n = 4, nev = 2, ncv = 4, ldv = n,
     &                  lworkl = ncv*(ncv+8))
      integer          iparam(11), ipntr(11), ido, info, ierr, i
      logical          select(ncv), rvec
      Real
     &                 resid(n), v(ldv,ncv), workd(3*n),
     &                 workl(lworkl), d(nev), z(ldv,nev), tol, sigma
c
      ido = 0
      info = 1
      tol = 1.0E-5
      sigma = 0.0E+0
      do 10 i = 1, n
         resid(i) = real(i)
   10 continue
      do 20 i = 1, 11
         iparam(i) = 0
   20 continue
      iparam(1) = 1
      iparam(3) = 100
      iparam(7) = 1
c
   30 continue
      call ssaupd ( ido, 'I', n, 'BE', nev, tol, resid, ncv, v, ldv,
     &              iparam, ipntr, workd, workl, lworkl, info )
      if (ido .eq. -1 .or. ido .eq. 1) then
         do 40 i = 1, n
            workd(ipntr(2)+i-1) = real(i) * workd(ipntr(1)+i-1)
   40    continue
         go to 30
      end if
c
      if (info .ne. 0 .or. iparam(5) .ne. nev) then
         print *, ' Error with _saupd, info = ', info
         stop 1
      end if
c
      rvec = .false.
      call sseupd ( rvec, 'A', select, d, z, ldv, sigma, 'I', n,
     &              'BE', nev, tol, resid, ncv, v, ldv, iparam,
     &              ipntr, workd, workl, lworkl, ierr )
      if (ierr .ne. 0) then
         print *, ' Error with _seupd, info = ', ierr
         stop 1
      end if
c
      print *, ' Ritz values: ', d(1), d(2)
      if (abs(d(1) - 1.0E+0) .gt. 1.0E-4 .or.
     &    abs(d(2) - 4.0E+0) .gt. 1.0E-4) then
         print *, ' Error: expected Ritz values 1 and 4'
         stop 1
      end if
      end
