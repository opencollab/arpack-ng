!\BeginDoc
!
!\Name: pznapps
!
! Message Passing Layer: MPI
!
!\Description:
!  Given the Arnoldi factorization
!
!     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!
!  apply NP implicit shifts resulting in
!
!     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!
!  where Q is an orthogonal matrix which is the product of rotations
!  and reflections resulting from the NP bulge change sweeps.
!  The updated Arnoldi factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call pznapps
!     ( COMM, N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ,
!       WORKL, WORKD )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  N       Integer.  (INPUT)
!          Problem size, i.e. size of matrix A.
!
!  KEV     Integer.  (INPUT/OUTPUT)
!          KEV+NP is the size of the input matrix H.
!          KEV is the size of the updated matrix HNEW.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFT   Complex*16 array of length NP.  (INPUT)
!          The shifts to be applied.
!
!  V       Complex*16 N by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, V contains the current KEV+NP Arnoldi vectors.
!          On OUTPUT, V contains the updated KEV Arnoldi vectors
!          in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Complex*16 (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, H contains the current KEV+NP by KEV+NP upper
!          Hessenberg matrix of the Arnoldi factorization.
!          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
!          matrix in the KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
!          On INPUT, RESID contains the the residual vector r_{k+p}.
!          On OUTPUT, RESID is the update residual vector rnew_{k}
!          in the first KEV locations.
!
!  Q       Complex*16 KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations and reflections
!          during the bulge chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Complex*16 work array of length (KEV+NP).  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  WORKD   Complex*16 work array of length 2*N.  (WORKSPACE)
!          Distributed array used in the application of the accumulated
!          orthogonal matrix Q.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  Complex*16
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!
!\Routines called:
!     pivout  Parallel ARPACK utility routine that prints integers.
!     arscnd  ARPACK utility routine for timing.
!     pzmout  Parallel ARPACK utility routine that prints matrices
!     pzvout  Parallel ARPACK utility routine that prints vectors.
!     zlacpy  LAPACK matrix copy routine.
!     zlanhs  LAPACK routine that computes various norms of a matrix.
!     zlartg  LAPACK Givens rotation construction routine.
!     zlaset  LAPACK matrix initialization routine.
!     dlabad  LAPACK routine for defining the underflow and overflow
!             limits.
!     pdlamch10 ScaLAPACK routine that determines machine constants.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     zgemv   Level 2 BLAS routine for matrix vector multiplication.
!     zaxpy   Level 1 BLAS that computes a vector triad.
!     zcopy   Level 1 BLAS that copies one vector to another.
!     zscal   Level 1 BLAS that scales a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Complex Code FILE: napps.F   SID: 2.1
!
!\SCCS Information:
! FILE: napps.F   SID: 1.4   DATE OF SID: 10/25/03
!
!\Remarks
!  1. In this version, each shift is applied to all the sublocks of
!     the Hessenberg matrix H and not just to the submatrix that it
!     comes from. Deflation as in LAPACK routine zlahqr (QR algorithm
!     for upper Hessenberg matrices ) is used.
!     Upon output, the subdiagonals of H are enforced to be non-negative
!     real numbers.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pznapps
     &   ( comm, n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq,
     &     workl, workd )
!
!     %--------------------%
!     | MPI Communicator |
!     %--------------------%
!
      include   'pcontext.h'
      integer   comm
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    kev, ldh, ldq, ldv, n, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Complex*16
     &           h(ldh,kev+np), resid(n), shift(np),
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex*16
     &           one, zero
      Double precision
     &           rzero
      parameter (one = (1.0, 0.0), zero = (0.0, 0.0),
     &           rzero = 0.0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      integer    i, iend, istart, j, jj, kplusp, msglvl
      Complex*16
     &           cdum, f, g, h11, h21, r, s, sigma, t
      Double precision
     &           c,  ovfl, smlnum, ulp, unfl, tst1
      save       ovfl, smlnum, ulp, unfl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   zaxpy, zcopy, zgemv, zscal, zlacpy, zlartg,
     &           pzvout, zlaset, dlabad, pzmout, arscnd, pivout
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &           zlanhs, pdlamch10, dlapy2
      external   zlanhs, pdlamch10, dlapy2
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs, dimag, conjg, dcmplx, max, min, dble
!
!     %---------------------%
!     | Statement Functions |
!     %---------------------%
!
      Double precision
     &           cabs1
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
!
!     %----------------%
!     | Data statements |
!     %----------------%
!
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (apps_first) then
!
!        %-----------------------------------------------%
!        | Set machine-dependent constants for the       |
!        | stopping criterion. If norm(H) <= sqrt(OVFL), |
!        | overflow should not occur.                    |
!        | REFERENCE: LAPACK subroutine zlahqr           |
!        %-----------------------------------------------%
!
         unfl = pdlamch10( comm, 'safe minimum' )
         ovfl = dble(one / unfl)
         call dlabad( unfl, ovfl )
         ulp = pdlamch10( comm, 'precision' )
         smlnum = unfl*( n / ulp )
         apps_first = .false.
      end if
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call arscnd (t0)
      msglvl = mcapps
!
      kplusp = kev + np
!
!     %--------------------------------------------%
!     | Initialize Q to the identity to accumulate |
!     | the rotations and reflections              |
!     %--------------------------------------------%
!
      call zlaset ('All', kplusp, kplusp, zero, one, q, ldq)
!
!     %----------------------------------------------%
!     | Quick return if there are no shifts to apply |
!     %----------------------------------------------%
!
      if (np .eq. 0) go to 9000
!
!     %----------------------------------------------%
!     | Chase the bulge with the application of each |
!     | implicit shift. Each shift is applied to the |
!     | whole matrix including each block.           |
!     %----------------------------------------------%
!
      do 110 jj = 1, np
         sigma = shift(jj)
!
         if (msglvl .gt. 2 ) then
            call pivout (comm, logfil, 1, [jj], ndigit,
     &               '_napps: shift number.')
            call pzvout (comm, logfil, 1, [sigma], ndigit,
     &               '_napps: Value of the shift ')
         end if
!
         istart = 1
   20    continue
!
         do 30 i = istart, kplusp-1
!
!           %----------------------------------------%
!           | Check for splitting and deflation. Use |
!           | a standard test as in the QR algorithm |
!           | REFERENCE: LAPACK subroutine zlahqr    |
!           %----------------------------------------%
!
            tst1 = cabs1( h( i, i ) ) + cabs1( h( i+1, i+1 ) )
            if( tst1.eq.rzero )
     &         tst1 = zlanhs( '1', kplusp-jj+1, h, ldh, workl )
            if ( abs(dble(h(i+1,i)))
     &           .le. max(ulp*tst1, smlnum) )  then
               if (msglvl .gt. 0) then
                  call pivout (comm, logfil, 1, [i], ndigit,
     &                 '_napps: matrix splitting at row/column no.')
                  call pivout (comm, logfil, 1, [jj], ndigit,
     &                 '_napps: matrix splitting with shift number.')
                  call pzvout (comm, logfil, 1, h(i+1,i), ndigit,
     &                 '_napps: off diagonal element.')
               end if
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
!
         if (msglvl .gt. 2) then
             call pivout (comm, logfil, 1, [istart], ndigit,
     &                   '_napps: Start of current block ')
             call pivout (comm, logfil, 1, [iend], ndigit,
     &                   '_napps: End of current block ')
         end if
!
!        %------------------------------------------------%
!        | No reason to apply a shift to block of order 1 |
!        | or if the current block starts after the point |
!        | of compression since we'll discard this stuff  |
!        %------------------------------------------------%
!
         if ( istart .eq. iend .or. istart .gt. kev) go to 100
!
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         f = h11 - sigma
         g = h21
!
         do 80 i = istart, iend-1
!
!           %------------------------------------------------------%
!           | Construct the plane rotation G to zero out the bulge |
!           %------------------------------------------------------%
!
            call zlartg (f, g, c, s, r)
            if (i .gt. istart) then
               h(i,i-1) = r
               h(i+1,i-1) = zero
            end if
!
!           %---------------------------------------------%
!           | Apply rotation to the left of H;  H <- G'*H |
!           %---------------------------------------------%
!
            do 50 j = i, kplusp
               t        =  c*h(i,j) + s*h(i+1,j)
               h(i+1,j) = -conjg(s)*h(i,j) + c*h(i+1,j)
               h(i,j)   = t
   50       continue
!
!           %---------------------------------------------%
!           | Apply rotation to the right of H;  H <- H*G |
!           %---------------------------------------------%
!
            do 60 j = 1, min(i+2,iend)
               t        =  c*h(j,i) + conjg(s)*h(j,i+1)
               h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
               h(j,i)   = t
   60       continue
!
!           %-----------------------------------------------------%
!           | Accumulate the rotation in the matrix Q;  Q <- Q*G' |
!           %-----------------------------------------------------%
!
            do 70 j = 1, min(i+jj, kplusp)
               t        =   c*q(j,i) + conjg(s)*q(j,i+1)
               q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
               q(j,i)   = t
   70       continue
!
!           %---------------------------%
!           | Prepare for next rotation |
!           %---------------------------%
!
            if (i .lt. iend-1) then
               f = h(i+1,i)
               g = h(i+2,i)
            end if
   80    continue
!
!        %-------------------------------%
!        | Finished applying the shift.  |
!        %-------------------------------%
!
  100    continue
!
!        %---------------------------------------------------------%
!        | Apply the same shift to the next block if there is any. |
!        %---------------------------------------------------------%
!
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
!
!        %---------------------------------------------%
!        | Loop back to the top to get the next shift. |
!        %---------------------------------------------%
!
  110 continue
!
!     %---------------------------------------------------%
!     | Perform a similarity transformation that makes    |
!     | sure that the compressed H will have non-negative |
!     | real subdiagonal elements.                        |
!     %---------------------------------------------------%
!
      do 120 j=1,kev
         if ( dble( h(j+1,j) ) .lt. rzero .or.
     &        dimag( h(j+1,j) ) .ne. rzero ) then
            t = h(j+1,j) / dlapy2(dble(h(j+1,j)),dimag(h(j+1,j)))
            call zscal( kplusp-j+1, conjg(t), h(j+1,j), ldh )
            call zscal( min(j+2, kplusp), t, h(1,j+1), 1 )
            call zscal( min(j+np+1,kplusp), t, q(1,j+1), 1 )
            h(j+1,j) = dcmplx( dble( h(j+1,j) ), rzero )
         end if
  120 continue
!
      do 130 i = 1, kev
!
!        %--------------------------------------------%
!        | Final check for splitting and deflation.   |
!        | Use a standard test as in the QR algorithm |
!        | REFERENCE: LAPACK subroutine zlahqr.       |
!        | Note: Since the subdiagonals of the        |
!        | compressed H are nonnegative real numbers, |
!        | we take advantage of this.                 |
!        %--------------------------------------------%
!
         tst1 = cabs1( h( i, i ) ) + cabs1( h( i+1, i+1 ) )
         if( tst1 .eq. rzero )
     &       tst1 = zlanhs( '1', kev, h, ldh, workl )
         if( dble( h( i+1,i ) ) .le. max( ulp*tst1, smlnum ) )
     &       h(i+1,i) = zero
 130  continue
!
!     %-------------------------------------------------%
!     | Compute the (kev+1)-st column of (V*Q) and      |
!     | temporarily store the result in WORKD(N+1:2*N). |
!     | This is needed in the residual update since we  |
!     | cannot GUARANTEE that the corresponding entry   |
!     | of H would be zero as in exact arithmetic.      |
!     %-------------------------------------------------%
!
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero,
     &                workd(n+1), 1)
!
!     %----------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order       |
!     | taking advantage of the upper Hessenberg structure of Q. |
!     %----------------------------------------------------------%
!
      do 140 i = 1, kev
         call zgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call zcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
      call zlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
!
!     %--------------------------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
!     %--------------------------------------------------------------%
!
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zcopy (n, workd(n+1), 1, v(1,kev+1), 1)
!
!     %-------------------------------------%
!     | Update the residual vector:         |
!     |    r <- sigmak*r + betak*v(:,kev+1) |
!     | where                               |
!     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
!     |    betak = e_{kev+1}'*H*e_{kev}     |
!     %-------------------------------------%
!
      call zscal (n, q(kplusp,kev), resid, 1)
      if ( dble( h(kev+1,kev) ) .gt. rzero )
     &   call zaxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
!
      if (msglvl .gt. 1) then
         call pzvout (comm, logfil, 1, q(kplusp,kev), ndigit,
     &        '_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}')
         call pzvout (comm, logfil, 1, h(kev+1,kev), ndigit,
     &        '_napps: betak = e_{kev+1}^T*H*e_{kev}')
         call pivout (comm, logfil, 1, [kev], ndigit,
     &               '_napps: Order of the final Hessenberg matrix ')
         if (msglvl .gt. 2) then
            call pzmout (comm, logfil, kev, kev, h, ldh, ndigit,
     &      '_napps: updated Hessenberg matrix H for next iteration')
         end if
!
      end if
!
 9000 continue
      call arscnd (t1)
      tcapps = tcapps + (t1 - t0)
!
      return
!
!     %----------------%
!     | End of pznapps |
!     %----------------%
!
      end
