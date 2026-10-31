      program default_tol
c
c     Nonpositive TOL is replaced by machine precision on the first
c     AUPD call and must be retained across reverse communication.
c
      implicit none
      call solve(.false., 0.0D+0)
      call solve(.false., -1.0D+0)
      call solve(.true., 0.0D+0)
      call solve(.true., -1.0D+0)
      end
c
      subroutine solve(upper, initial_tol)
      implicit none
      logical          upper, select(6)
      integer          n, nev, ncv, lworkl
      parameter        (n = 10, nev = 2, ncv = 6,
     &                  lworkl = 3*ncv*(ncv+2))
      integer          ido, info, iparam(11), ipntr(14), i, j
      double precision initial_tol, tol, eps, dlamch, expected,
     &                 resid(n), v(n,ncv), workd(3*n), workl(lworkl),
     &                 workev(3*ncv), dr(nev+1), di(nev+1),
     &                 z(n,nev+1), residual2, norm2, r
      external         dlamch
c
c     A has diagonal 1,...,10 and, if UPPER, superdiagonal 1.
c     Use a deterministic starting vector for both solvers.
c
      ido = 0
      info = 1
      tol = initial_tol
      eps = dlamch('EpsMach')
      do 10 i = 1, n
         resid(i) = dble(i)
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
         call dnaupd(ido, 'I', n, 'LM', nev, tol, resid, ncv, v, n,
     &               iparam, ipntr, workd, workl, lworkl, info)
      else
         call dsaupd(ido, 'I', n, 'LM', nev, tol, resid, ncv, v, n,
     &               iparam, ipntr, workd, workl, lworkl, info)
      end if
      if (tol .ne. eps) then
         print *, ' AUPD did not retain default TOL: ', tol
         stop 1
      end if
      if (ido .eq. -1 .or. ido .eq. 1) then
         do 40 i = 1, n
            workd(ipntr(2)+i-1) = dble(i)*workd(ipntr(1)+i-1)
            if (upper .and. i .lt. n) then
               workd(ipntr(2)+i-1) = workd(ipntr(2)+i-1)
     &                             + workd(ipntr(1)+i)
            end if
   40    continue
         go to 30
      end if
      if (ido .ne. 99 .or. info .ne. 0 .or.
     &    iparam(5) .ne. nev) then
         print *, ' AUPD: upper, initial TOL, IDO, INFO, NCONV = ',
     &            upper, initial_tol, ido, info, iparam(5)
         stop 1
      end if
c
c     Pass the same TOL variable to EUPD, as in the Fortran examples.
c
      if (upper) then
         call dneupd(.true., 'A', select, dr, di, z, n, 0.0D+0,
     &               0.0D+0, workev, 'I', n, 'LM', nev, tol, resid,
     &               ncv, v, n, iparam, ipntr, workd, workl,
     &               lworkl, info)
      else
         call dseupd(.true., 'A', select, dr, z, n, 0.0D+0, 'I', n,
     &               'LM', nev, tol, resid, ncv, v, n, iparam, ipntr,
     &               workd, workl, lworkl, info)
         di(1) = 0.0D+0
         di(2) = 0.0D+0
      end if
      if (info .ne. 0 .or. tol .ne. eps) then
         print *, ' EUPD: INFO, TOL = ', info, tol
         stop 1
      end if
c
c     Check the two largest eigenvalues, vector norms and residuals
c     with the same bounds as icb_default_tol.c.
c
      do 60 j = 1, nev
         expected = dble(n-nev+j)
         if (upper .and. dr(1) .gt. dr(2)) expected = dble(n+1-j)
         residual2 = 0.0D+0
         norm2 = 0.0D+0
         do 50 i = 1, n
            r = (dble(i)-dr(j))*z(i,j)
            if (upper .and. i .lt. n) r = r + z(i+1,j)
            residual2 = residual2 + r*r
            norm2 = norm2 + z(i,j)*z(i,j)
   50    continue
         if (.not. (abs(dr(j)-expected) .le. 1.0D-10 .and.
     &              abs(di(j)) .le. 1.0D-10 .and.
     &              norm2 .ge. 0.9D+0 .and.
     &              residual2 .le. 1.0D-20*norm2)) then
            print *, ' Incorrect eigenpair: ', j, dr(j), di(j),
     &               norm2, residual2
            stop 1
         end if
   60 continue
      print *, ' upper, initial TOL, Ritz values: ', upper,
     &         initial_tol, dr(1), dr(2)
      end
