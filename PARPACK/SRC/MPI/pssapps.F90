!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: pssapps
!
!\Description:
!  Given the Arnoldi factorization
!
!     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!
!  apply NP shifts implicitly resulting in
!
!     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!
!  where Q is an orthogonal matrix of order KEV+NP. Q is the product of
!  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi
!  factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call pssapps
!     ( COMM, N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  N       Integer.  (INPUT)
!          Problem size, i.e. dimension of matrix A.
!
!  KEV     Integer.  (INPUT)
!          INPUT: KEV+NP is the size of the input matrix H.
!          OUTPUT: KEV is the size of the updated matrix HNEW.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFT   Real array of length NP.  (INPUT)
!          The shifts to be applied.
!
!  V       Real N by (KEV+NP) array.  (INPUT/OUTPUT)
!          INPUT: V contains the current KEV+NP Arnoldi vectors.
!          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
!          are in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (KEV+NP) by 2 array.  (INPUT/OUTPUT)
!          INPUT: H contains the symmetric tridiagonal matrix of the
!          Arnoldi factorization with the subdiagonal in the 1st column
!          starting at H(2,1) and the main diagonal in the 2nd column.
!          OUTPUT: H contains the updated tridiagonal matrix in the
!          KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   Real array of length (N).  (INPUT/OUTPUT)
!          INPUT: RESID contains the the residual vector r_{k+p}.
!          OUTPUT: RESID is the updated residual vector rnew_{k}.
!
!  Q       Real KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations during the bulge
!          chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
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
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     pivout  Parallel ARPACK utility routine that prints integers.
!     arscnd  ARPACK utility routine for timing.
!     psvout  Parallel ARPACK utility routine that prints vectors.
!     pslamch10 ScaLAPACK routine that determines machine constants.
!     slartg  LAPACK Givens rotation construction routine.
!     slacpy  LAPACK matrix copy routine.
!     slaset  LAPACK matrix initialization routine.
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     saxpy   Level 1 BLAS that computes a vector triad.
!     scopy   Level 1 BLAS that copies one vector to another.
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
!     Starting Point: Serial Code FILE: sapps.F   SID: 2.4
!
!\SCCS Information:
! FILE: sapps.F   SID: 1.3   DATE OF SID: 3/19/97
!
!\Remarks
!  1. In this version, each shift is applied to all the subblocks of
!     the tridiagonal matrix H and not just to the submatrix that it
!     comes from. This routine assumes that the subdiagonal elements
!     of H that are stored in h(1:kev+np,1) are nonegative upon input
!     and enforce this condition upon output. This version incorporates
!     deflation. See code for documentation.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pssapps
     &  ( comm, n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd)
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
     &           h(ldh,2), q(ldq,kev+np), resid(n), shift(np),
     &           v(ldv,kev+np), workd(2*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, iend, istart, itop, j, jj, kplusp, msglvl
      Real
     &           a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      save       epsmch
!
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   saxpy, scopy, sscal, slacpy, slartg, slaset, psvout,
     &           pivout, arscnd, sgemv
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           pslamch10
      external   pslamch10
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs
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
         epsmch = pslamch10(comm, 'Epsilon-Machine')
         apps_first = .false.
      end if
      itop = 1
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call arscnd (t0)
      msglvl = msapps
!
      kplusp = kev + np
!
!     %----------------------------------------------%
!     | Initialize Q to the identity matrix of order |
!     | kplusp used to accumulate the rotations.     |
!     %----------------------------------------------%
!
      call slaset ('All', kplusp, kplusp, zero, one, q, ldq)
!
!     %----------------------------------------------%
!     | Quick return if there are no shifts to apply |
!     %----------------------------------------------%
!
      if (np .eq. 0) go to 9000
!
!     %----------------------------------------------------------%
!     | Apply the np shifts implicitly. Apply each shift to the  |
!     | whole matrix and not just to the submatrix from which it |
!     | comes.                                                   |
!     %----------------------------------------------------------%
!
      do 90 jj = 1, np
!
         istart = itop
!
!        %----------------------------------------------------------%
!        | Check for splitting and deflation. Currently we consider |
!        | an off-diagonal element h(i+1,1) negligible if           |
!        |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
!        | for i=1:KEV+NP-1.                                        |
!        | If above condition tests true then we set h(i+1,1) = 0.  |
!        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
!        %----------------------------------------------------------%
!
   20    continue
!
!        %------------------------------------------------%
!        | The following loop exits early if we encounter |
!        | a negligible off diagonal element.             |
!        %------------------------------------------------%
!
         do 30 i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) .le. epsmch*big) then
               if (msglvl .gt. 0) then
                  call pivout (comm, logfil, 1, [i], ndigit,
     &                 '_sapps: deflation at row/column no.')
                  call pivout (comm, logfil, 1, [jj], ndigit,
     &                 '_sapps: occurred before shift number.')
                  call psvout (comm, logfil, 1, h(i+1,1), ndigit,
     &                 '_sapps: the corresponding off diagonal element')
               end if
               h(i+1,1) = zero
               iend = i
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
!
         if (istart .lt. iend) then
!
!           %--------------------------------------------------------%
!           | Construct the plane rotation G'(istart,istart+1,theta) |
!           | that attempts to drive h(istart+1,1) to zero.          |
!           %--------------------------------------------------------%
!
             f = h(istart,2) - shift(jj)
             g = h(istart+1,1)
             call slartg (f, g, c, s, r)
!
!            %-------------------------------------------------------%
!            | Apply rotation to the left and right of H;            |
!            | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
!            | This will create a "bulge".                           |
!            %-------------------------------------------------------%
!
             a1 = c*h(istart,2)   + s*h(istart+1,1)
             a2 = c*h(istart+1,1) + s*h(istart+1,2)
             a4 = c*h(istart+1,2) - s*h(istart+1,1)
             a3 = c*h(istart+1,1) - s*h(istart,2)
             h(istart,2)   = c*a1 + s*a2
             h(istart+1,2) = c*a4 - s*a3
             h(istart+1,1) = c*a3 + s*a4
!
!            %----------------------------------------------------%
!            | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!            %----------------------------------------------------%
!
             do 60 j = 1, min(istart+jj,kplusp)
                a1            =   c*q(j,istart) + s*q(j,istart+1)
                q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
                q(j,istart)   = a1
   60        continue
!
!
!            %----------------------------------------------%
!            | The following loop chases the bulge created. |
!            | Note that the previous rotation may also be  |
!            | done within the following loop. But it is    |
!            | kept separate to make the distinction among  |
!            | the bulge chasing sweeps and the first plane |
!            | rotation designed to drive h(istart+1,1) to  |
!            | zero.                                        |
!            %----------------------------------------------%
!
             do 70 i = istart+1, iend-1
!
!               %----------------------------------------------%
!               | Construct the plane rotation G'(i,i+1,theta) |
!               | that zeros the i-th bulge that was created   |
!               | by G(i-1,i,theta). g represents the bulge.   |
!               %----------------------------------------------%
!
                f = h(i,1)
                g = s*h(i+1,1)
!
!               %----------------------------------%
!               | Final update with G(i-1,i,theta) |
!               %----------------------------------%
!
                h(i+1,1) = c*h(i+1,1)
                call slartg (f, g, c, s, r)
!
!               %-------------------------------------------%
!               | The following ensures that h(1:iend-1,1), |
!               | the first iend-2 off diagonal of elements |
!               | H, remain non negative.                   |
!               %-------------------------------------------%
!
                if (r .lt. zero) then
                   r = -r
                   c = -c
                   s = -s
                end if
!
!               %--------------------------------------------%
!               | Apply rotation to the left and right of H; |
!               | H <- G * H * G',  where G = G(i,i+1,theta) |
!               %--------------------------------------------%
!
                h(i,1) = r
!
                a1 = c*h(i,2)   + s*h(i+1,1)
                a2 = c*h(i+1,1) + s*h(i+1,2)
                a3 = c*h(i+1,1) - s*h(i,2)
                a4 = c*h(i+1,2) - s*h(i+1,1)
!
                h(i,2)   = c*a1 + s*a2
                h(i+1,2) = c*a4 - s*a3
                h(i+1,1) = c*a3 + s*a4
!
!               %----------------------------------------------------%
!               | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!               %----------------------------------------------------%
!
                do 50 j = 1, min( i+jj, kplusp )
                   a1       =   c*q(j,i) + s*q(j,i+1)
                   q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                   q(j,i)   = a1
   50           continue
!
   70        continue
!
         end if
!
!        %--------------------------%
!        | Update the block pointer |
!        %--------------------------%
!
         istart = iend + 1
!
!        %------------------------------------------%
!        | Make sure that h(iend,1) is non-negative |
!        | If not then set h(iend,1) <-- -h(iend,1) |
!        | and negate the last column of Q.         |
!        | We have effectively carried out a        |
!        | similarity on transformation H           |
!        %------------------------------------------%
!
         if (h(iend,1) .lt. zero) then
             h(iend,1) = -h(iend,1)
             call sscal(kplusp, -one, q(1,iend), 1)
         end if
!
!        %--------------------------------------------------------%
!        | Apply the same shift to the next block if there is any |
!        %--------------------------------------------------------%
!
         if (iend .lt. kplusp) go to 20
!
!        %-----------------------------------------------------%
!        | Check if we can increase the the start of the block |
!        %-----------------------------------------------------%
!
         do 80 i = itop, kplusp-1
            if (h(i+1,1) .gt. zero) go to 90
            itop  = itop + 1
   80    continue
!
!        %-----------------------------------%
!        | Finished applying the jj-th shift |
!        %-----------------------------------%
!
   90 continue
!
!     %------------------------------------------%
!     | All shifts have been applied. Check for  |
!     | more possible deflation that might occur |
!     | after the last shift is applied.         |
!     %------------------------------------------%
!
      do 100 i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) .le. epsmch*big) then
            if (msglvl .gt. 0) then
               call pivout (comm, logfil, 1, [i], ndigit,
     &              '_sapps: deflation at row/column no.')
               call psvout (comm, logfil, 1, h(i+1,1), ndigit,
     &              '_sapps: the corresponding off diagonal element')
            end if
            h(i+1,1) = zero
         end if
 100  continue
!
!     %-------------------------------------------------%
!     | Compute the (kev+1)-st column of (V*Q) and      |
!     | temporarily store the result in WORKD(N+1:2*N). |
!     | This is not necessary if h(kev+1,1) = 0.         |
!     %-------------------------------------------------%
!
      if ( h(kev+1,1) .gt. zero )
     &   call sgemv ('N', n, kplusp, one, v, ldv,
     &                q(1,kev+1), 1, zero, workd(n+1), 1)
!
!     %-------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order    |
!     | taking advantage that Q is an upper triangular matrix |
!     | with lower bandwidth np.                              |
!     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
!     %-------------------------------------------------------%
!
      do 130 i = 1, kev
         call sgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call scopy (n, workd, 1, v(1,kplusp-i+1), 1)
  130 continue
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
      call slacpy ('All', n, kev, v(1,np+1), ldv, v, ldv)
!
!     %--------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the |
!     | appropriate place if h(kev+1,1) .ne. zero. |
!     %--------------------------------------------%
!
      if ( h(kev+1,1) .gt. zero )
     &     call scopy (n, workd(n+1), 1, v(1,kev+1), 1)
!
!     %-------------------------------------%
!     | Update the residual vector:         |
!     |    r <- sigmak*r + betak*v(:,kev+1) |
!     | where                               |
!     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
!     |    betak = e_{kev+1}'*H*e_{kev}     |
!     %-------------------------------------%
!
      call sscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,1) .gt. zero)
     &   call saxpy (n, h(kev+1,1), v(1,kev+1), 1, resid, 1)
!
      if (msglvl .gt. 1) then
         call psvout (comm, logfil, 1, q(kplusp,kev), ndigit,
     &      '_sapps: sigmak of the updated residual vector')
         call psvout (comm, logfil, 1, h(kev+1,1), ndigit,
     &      '_sapps: betak of the updated residual vector')
         call psvout (comm, logfil, kev, h(1,2), ndigit,
     &      '_sapps: updated main diagonal of H for next iteration')
         if (kev .gt. 1) then
         call psvout (comm, logfil, kev-1, h(2,1), ndigit,
     &      '_sapps: updated sub diagonal of H for next iteration')
         end if
      end if
!
      call arscnd (t1)
      tsapps = tsapps + (t1 - t0)
!
 9000 continue
      return
!
!     %----------------%
!     | End of pssapps |
!     %----------------%
!
      end
