!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: psnapps
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
!  and reflections resulting from the NP bulge chage sweeps.
!  The updated Arnoldi factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call psnapps
!     ( COMM, N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ,
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
!          KEV is the size of the updated matrix HNEW.  KEV is only
!          updated on output when fewer than NP shifts are applied in
!          order to keep the conjugate pair together.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFTR, Real array of length NP.  (INPUT)
!  SHIFTI  Real and imaginary part of the shifts to be applied.
!          Upon, entry to psnapps, the shifts must be sorted so that the
!          conjugate pairs are in consecutive locations.
!
!  V       Real N by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, V contains the current KEV+NP Arnoldi vectors.
!          On OUTPUT, V contains the updated KEV Arnoldi vectors
!          in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, H contains the current KEV+NP by KEV+NP upper
!          Hessenber matrix of the Arnoldi factorization.
!          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
!          matrix in the KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          On INPUT, RESID contains the the residual vector r_{k+p}.
!          On OUTPUT, RESID is the update residual vector rnew_{k}
!          in the first KEV locations.
!
!  Q       Real KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations and reflections
!          during the bulge chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Real work array of length (KEV+NP).  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  WORKD   Real work array of length 2*N.  (WORKSPACE)
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
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!
!\Routines called:
!     slabad  LAPACK routine that computes machine constants.
!     slacpy  LAPACK matrix copy routine.
!     pslamch10 ScaLAPACK routine that determines machine constants.
!     slanhs  LAPACK routine that computes various norms of a matrix.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     slarf   LAPACK routine that applies Householder reflection to
!             a matrix.
!     slarfg  LAPACK Householder reflection construction routine.
!     slartg  LAPACK Givens rotation construction routine.
!     slaset  LAPACK matrix initialization routine.
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     saxpy   Level 1 BLAS that computes a vector triad.
!     scopy   Level 1 BLAS that copies one vector to another .
!     sscal   Level 1 BLAS that scales a vector.
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
!     Starting Point: Serial Code FILE: napps.F   SID: 2.2
!
!\SCCS Information:
! FILE: napps.F   SID: 1.5   DATE OF SID: 03/19/97
!
!\Remarks
!  1. In this version, each shift is applied to all the sublocks of
!     the Hessenberg matrix H and not just to the submatrix that it
!     comes from. Deflation as in LAPACK routine slahqr (QR algorithm
!     for upper Hessenberg matrices ) is used.
!     The subdiagonals of H are enforced to be non-negative.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine psnapps
     &   ( comm, n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid,
     &     q, ldq, workl, workd )
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
      Real
     &           h(ldh,kev+np), resid(n), shifti(np), shiftr(np),
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      integer    i, iend, ir, istart, j, jj, kplusp, msglvl, nr
      logical    cconj
      Real
     &           c, f, g, h11, h12, h21, h22, h32, ovfl, r, s, sigmai,
     &           sigmar, smlnum, ulp, unfl, u(3), t, tau, tst1
      save       ovfl, smlnum, ulp, unfl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   saxpy, scopy, sscal, slacpy, slarf, slarfg, slartg,
     &           slaset, slabad, arscnd, pivout, psvout, psmout
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           pslamch10, slanhs, slapy2
      external   pslamch10, slanhs, slapy2
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs, max, min
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
!        | REFERENCE: LAPACK subroutine slahqr           |
!        %-----------------------------------------------%
!
         unfl = pslamch10( comm, 'safe minimum' )
         ovfl = one / unfl
         call slabad( unfl, ovfl )
         ulp = pslamch10( comm, 'precision' )
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
      msglvl = mnapps
!
      kplusp = kev + np
!
!     %--------------------------------------------%
!     | Initialize Q to the identity to accumulate |
!     | the rotations and reflections              |
!     %--------------------------------------------%
!
      call slaset ('All', kplusp, kplusp, zero, one, q, ldq)
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
      cconj = .false.
      do 110 jj = 1, np
         sigmar = shiftr(jj)
         sigmai = shifti(jj)
!
         if (msglvl .gt. 2 ) then
            call pivout (comm, logfil, 1, [jj], ndigit,
     &               '_napps: shift number.')
            call psvout (comm, logfil, 1, [sigmar], ndigit,
     &               '_napps: The real part of the shift ')
            call psvout (comm, logfil, 1, [sigmai], ndigit,
     &               '_napps: The imaginary part of the shift ')
         end if
!
!        %-------------------------------------------------%
!        | The following set of conditionals is necessary  |
!        | in order that complex conjugate pairs of shifts |
!        | are applied together or not at all.             |
!        %-------------------------------------------------%
!
         if ( cconj ) then
!
!           %-----------------------------------------%
!           | cconj = .true. means the previous shift |
!           | had non-zero imaginary part.            |
!           %-----------------------------------------%
!
            cconj = .false.
            go to 110
         else if ( jj .lt. np .and. abs( sigmai ) .gt. zero ) then
!
!           %------------------------------------%
!           | Start of a complex conjugate pair. |
!           %------------------------------------%
!
            cconj = .true.
         else if ( jj .eq. np .and. abs( sigmai ) .gt. zero ) then
!
!           %----------------------------------------------%
!           | The last shift has a nonzero imaginary part. |
!           | Don't apply it; thus the order of the        |
!           | compressed H is order KEV+1 since only np-1  |
!           | were applied.                                |
!           %----------------------------------------------%
!
            kev = kev + 1
            go to 110
         end if
         istart = 1
   20    continue
!
!        %--------------------------------------------------%
!        | if sigmai = 0 then                               |
!        |    Apply the jj-th shift ...                     |
!        | else                                             |
!        |    Apply the jj-th and (jj+1)-th together ...    |
!        |    (Note that jj < np at this point in the code) |
!        | end                                              |
!        | to the current block of H. The next do loop      |
!        | determines the current block ;                   |
!        %--------------------------------------------------%
!
         do 30 i = istart, kplusp-1
!
!           %----------------------------------------%
!           | Check for splitting and deflation. Use |
!           | a standard test as in the QR algorithm |
!           | REFERENCE: LAPACK subroutine slahqr    |
!           %----------------------------------------%
!
            tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
            if( tst1.eq.zero )
     &         tst1 = slanhs( '1', kplusp-jj+1, h, ldh, workl )
            if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) then
               if (msglvl .gt. 0) then
                  call pivout (comm, logfil, 1, [i], ndigit,
     &                 '_napps: matrix splitting at row/column no.')
                  call pivout (comm, logfil, 1, [jj], ndigit,
     &                 '_napps: matrix splitting with shift number.')
                  call psvout (comm, logfil, 1, h(i+1,i), ndigit,
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
!        %------------------------------------------------%
!
         if ( istart .eq. iend ) go to 100
!
!        %------------------------------------------------------%
!        | If istart + 1 = iend then no reason to apply a       |
!        | complex conjugate pair of shifts on a 2 by 2 matrix. |
!        %------------------------------------------------------%
!
         if ( istart + 1 .eq. iend .and. abs( sigmai ) .gt. zero )
     &      go to 100
!
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         if ( abs( sigmai ) .le. zero ) then
!
!           %---------------------------------------------%
!           | Real-valued shift ==> apply single shift QR |
!           %---------------------------------------------%
!
            f = h11 - sigmar
            g = h21
!
            do 80 i = istart, iend-1
!
!              %-----------------------------------------------------%
!              | Construct the plane rotation G to zero out the bulge |
!              %-----------------------------------------------------%
!
               call slartg (f, g, c, s, r)
               if (i .gt. istart) then
!
!                 %-------------------------------------------%
!                 | The following ensures that h(1:iend-1,1), |
!                 | the first iend-2 off diagonal of elements |
!                 | H, remain non negative.                   |
!                 %-------------------------------------------%
!
                  if (r .lt. zero) then
                     r = -r
                     c = -c
                     s = -s
                  end if
                  h(i,i-1) = r
                  h(i+1,i-1) = zero
               end if
!
!              %---------------------------------------------%
!              | Apply rotation to the left of H;  H <- G'*H |
!              %---------------------------------------------%
!
               do 50 j = i, kplusp
                  t        =  c*h(i,j) + s*h(i+1,j)
                  h(i+1,j) = -s*h(i,j) + c*h(i+1,j)
                  h(i,j)   = t
   50          continue
!
!              %---------------------------------------------%
!              | Apply rotation to the right of H;  H <- H*G |
!              %---------------------------------------------%
!
               do 60 j = 1, min(i+2,iend)
                  t        =  c*h(j,i) + s*h(j,i+1)
                  h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
                  h(j,i)   = t
   60          continue
!
!              %----------------------------------------------------%
!              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!              %----------------------------------------------------%
!
               do 70 j = 1, min( i+jj, kplusp )
                  t        =   c*q(j,i) + s*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = t
   70          continue
!
!              %---------------------------%
!              | Prepare for next rotation |
!              %---------------------------%
!
               if (i .lt. iend-1) then
                  f = h(i+1,i)
                  g = h(i+2,i)
               end if
   80       continue
!
!           %-----------------------------------%
!           | Finished applying the real shift. |
!           %-----------------------------------%
!
         else
!
!           %----------------------------------------------------%
!           | Complex conjugate shifts ==> apply double shift QR |
!           %----------------------------------------------------%
!
            h12 = h(istart,istart+1)
            h22 = h(istart+1,istart+1)
            h32 = h(istart+2,istart+1)
!
!           %---------------------------------------------------------%
!           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
!           %---------------------------------------------------------%
!
            s    = 2.0*sigmar
            t = slapy2 ( sigmar, sigmai )
            u(1) = ( h11 * (h11 - s) + t * t ) / h21 + h12
            u(2) = h11 + h22 - s
            u(3) = h32
!
            do 90 i = istart, iend-1
!
               nr = min ( 3, iend-i+1 )
!
!              %-----------------------------------------------------%
!              | Construct Householder reflector G to zero out u(1). |
!              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
!              %-----------------------------------------------------%
!
               call slarfg ( nr, u(1), u(2), 1, tau )
!
               if (i .gt. istart) then
                  h(i,i-1)   = u(1)
                  h(i+1,i-1) = zero
                  if (i .lt. iend-1) h(i+2,i-1) = zero
               end if
               u(1) = one
!
!              %--------------------------------------%
!              | Apply the reflector to the left of H |
!              %--------------------------------------%
!
               call slarf ('Left', nr, kplusp-i+1, u, 1, tau,
     &                     h(i,i), ldh, workl)
!
!              %---------------------------------------%
!              | Apply the reflector to the right of H |
!              %---------------------------------------%
!
               ir = min ( i+3, iend )
               call slarf ('Right', ir, nr, u, 1, tau,
     &                     h(1,i), ldh, workl)
!
!              %-----------------------------------------------------%
!              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
!              %-----------------------------------------------------%
!
               call slarf ('Right', kplusp, nr, u, 1, tau,
     &                     q(1,i), ldq, workl)
!
!              %----------------------------%
!              | Prepare for next reflector |
!              %----------------------------%
!
               if (i .lt. iend-1) then
                  u(1) = h(i+1,i)
                  u(2) = h(i+2,i)
                  if (i .lt. iend-2) u(3) = h(i+3,i)
               end if
!
   90       continue
!
!           %--------------------------------------------%
!           | Finished applying a complex pair of shifts |
!           | to the current block                       |
!           %--------------------------------------------%
!
         end if
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
!     %--------------------------------------------------%
!     | Perform a similarity transformation that makes   |
!     | sure that H will have non negative sub diagonals |
!     %--------------------------------------------------%
!
      do 120 j=1,kev
         if ( h(j+1,j) .lt. zero ) then
              call sscal( kplusp-j+1, -one, h(j+1,j), ldh )
              call sscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
              call sscal( min(j+np+1,kplusp), -one, q(1,j+1), 1 )
         end if
 120  continue
!
      do 130 i = 1, kev
!
!        %--------------------------------------------%
!        | Final check for splitting and deflation.   |
!        | Use a standard test as in the QR algorithm |
!        | REFERENCE: LAPACK subroutine slahqr        |
!        %--------------------------------------------%
!
         tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
         if( tst1.eq.zero )
     &       tst1 = slanhs( '1', kev, h, ldh, workl )
         if( h( i+1,i ) .le. max( ulp*tst1, smlnum ) )
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
      if (h(kev+1,kev) .gt. zero)
     &    call sgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero,
     &                workd(n+1), 1)
!
!     %----------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order       |
!     | taking advantage of the upper Hessenberg structure of Q. |
!     %----------------------------------------------------------%
!
      do 140 i = 1, kev
         call sgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call scopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
      call slacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
!
!     %--------------------------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
!     %--------------------------------------------------------------%
!
      if (h(kev+1,kev) .gt. zero)
     &   call scopy (n, workd(n+1), 1, v(1,kev+1), 1)
!
!     %---------------------------------------%
!     | Update the residual vector:           |
!     |    r <- sigmak*r + betak*v(:,kev+1)   |
!     | where                                 |
!     |    sigmak = (e_{kplusp}'*Q)*e_{kev}   |
!     |    betak = e_{kev+1}'*H*e_{kev}       |
!     %---------------------------------------%
!
      call sscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,kev) .gt. zero)
     &   call saxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
!
      if (msglvl .gt. 1) then
         call psvout (comm, logfil, 1, q(kplusp,kev), ndigit,
     &        '_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}')
         call psvout (comm, logfil, 1, h(kev+1,kev), ndigit,
     &        '_napps: betak = e_{kev+1}^T*H*e_{kev}')
         call pivout (comm, logfil, 1, [kev], ndigit,
     &               '_napps: Order of the final Hessenberg matrix ')
         if (msglvl .gt. 2) then
            call psmout (comm, logfil, kev, kev, h, ldh, ndigit,
     &      '_napps: updated Hessenberg matrix H for next iteration')
         end if
!
      end if
!
 9000 continue
      call arscnd (t1)
      tnapps = tnapps + (t1 - t0)
!
      return
!
!     %----------------%
!     | End of psnapps |
!     %----------------%
!
      end
