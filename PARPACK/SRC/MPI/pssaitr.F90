!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: pssaitr
!
! Message Passing Layer: MPI
!
!\Description:
!  Reverse communication interface for applying NP additional steps to
!  a K step symmetric Arnoldi factorization.
!
!  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
!
!          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
!
!  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
!
!          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
!
!  where OP and B are as in pssaupd.  The B-norm of r_{k+p} is also
!  computed and returned.
!
!\Usage:
!  call pssaitr
!     ( COMM, IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH,
!       IPNTR, WORKD, WORKL, INFO )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!                    This is for the restart phase to force the new
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y,
!                    IPNTR(3) is the pointer into WORK for B * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!          When the routine is used in the "shift-and-invert" mode, the
!          vector B * Q is already available and does not need to be
!          recomputed in forming OP * Q.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of matrix B that defines the
!          semi-inner product for the operator OP.  See pssaupd.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  K       Integer.  (INPUT)
!          Current order of H and the number of columns of V.
!
!  NP      Integer.  (INPUT)
!          Number of additional Arnoldi steps to take.
!
!  MODE    Integer.  (INPUT)
!          Signifies which form for "OP". If MODE=2 then
!          a reduction in the number of B matrix vector multiplies
!          is possible since the B-norm of OP*x is equivalent to
!          the inv(B)-norm of A*x.
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          On INPUT:  RESID contains the residual vector r_{k}.
!          On OUTPUT: RESID contains the residual vector r_{k+p}.
!
!  RNORM   Real scalar.  (INPUT/OUTPUT)
!          On INPUT the B-norm of r_{k}.
!          On OUTPUT the B-norm of the updated residual r_{k+p}.
!
!  V       Real N by K+NP array.  (INPUT/OUTPUT)
!          On INPUT:  V contains the Arnoldi vectors in the first K
!          columns.
!          On OUTPUT: V contains the new NP Arnoldi vectors in the next
!          NP columns.  The first K columns are unchanged.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (K+NP) by 2 array.  (INPUT/OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          with the subdiagonal in the first column starting at H(2,1)
!          and the main diagonal in the second column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORK for
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!
!  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The calling program should not
!          use WORKD as temporary workspace during the iteration !!!!!!
!          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K step Arnoldi factorization. Used to save some
!          computation at the first step.
!          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K+NP step Arnoldi factorization.
!
!  WORKL   Real work space used for Gram Schmidt orthogonalization
!
!  INFO    Integer.  (OUTPUT)
!          = 0: Normal exit.
!          > 0: Size of an invariant subspace of OP is found that is
!               less than K + NP.
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
!\Routines called:
!     psgetv0  Parallel ARPACK routine to generate the initial vector.
!     pivout   Parallel ARPACK utility routine that prints integers.
!     psmout   Parallel ARPACK utility routine that prints matrices.
!     psvout   Parallel ARPACK utility routine that prints vectors.
!     pslamch10  ScaLAPACK routine that determines machine constants.
!     slascl   LAPACK routine for careful scaling of a matrix.
!     sgemv    Level 2 BLAS routine for matrix vector multiplication.
!     saxpy    Level 1 BLAS that computes a vector triad.
!     sscal    Level 1 BLAS that scales a vector.
!     scopy    Level 1 BLAS that copies one vector to another .
!     sdot     Level 1 BLAS that computes the scalar product of two vectors.
!     psnorm2  Parallel version of Level 1 BLAS that computes the norm of a vector.
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
!     Starting Point: Serial Code FILE: saitr.F   SID: 2.3
!
!\SCCS Information:
! FILE: saitr.F   SID: 1.3   DATE OF SID: 3/19/97
!
!\Remarks
!  The algorithm implemented is:
!
!  restart = .false.
!  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
!  r_{k} contains the initial residual vector even for k = 0;
!  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
!  computed by the calling program.
!
!  betaj = rnorm ; p_{k+1} = B*r_{k} ;
!  For  j = k+1, ..., k+np  Do
!     1) if ( betaj < tol ) stop or restart depending on j.
!        if ( restart ) generate a new starting vector.
!     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
!        p_{j} = p_{j}/betaj
!     3) r_{j} = OP*v_{j} where OP is defined as in pssaupd
!        For shift-invert mode p_{j} = B*v_{j} is already available.
!        wnorm = || OP*v_{j} ||
!     4) Compute the j-th step residual vector.
!        w_{j} =  V_{j}^T * B * OP * v_{j}
!        r_{j} =  OP*v_{j} - V_{j} * w_{j}
!        alphaj <- j-th component of w_{j}
!        rnorm = || r_{j} ||
!        betaj+1 = rnorm
!        If (rnorm > 0.717*wnorm) accept step and go back to 1)
!     5) Re-orthogonalization step:
!        s = V_{j}'*B*r_{j}
!        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
!        alphaj = alphaj + s_{j};
!     6) Iterative refinement step:
!        If (rnorm1 > 0.717*rnorm) then
!           rnorm = rnorm1
!           accept step and go back to 1)
!        Else
!           rnorm = rnorm1
!           If this is the first time in step 6), go to 5)
!           Else r_{j} lies in the span of V_{j} numerically.
!              Set r_{j} = 0 and rnorm = 0; go to 1)
!        EndIf
!  End Do
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine pssaitr
     &   (comm, ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh,
     &    ipntr, workd, workl, info)
!
      include   'pcontext.h'
      include   'mpif.h'
!
!     %---------------%
!     | MPI Variables |
!     %---------------%
!
      integer    comm
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
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, mode, np
      Real
     &           rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Real
     &           h(ldh,2), resid(n), v(ldv,k+np), workd(3*n),
     &           workl(2*ldh)
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
      logical    orth1, orth2, rstart, step3, step4
      integer    i, ierr, ipj, irj, ivj, iter, itry, j, msglvl, infol,
     &           jj
      Real
     &           rnorm1, wnorm(1), safmin, temp1, temp2(1)
      save       orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl,
     &           rnorm1, safmin, wnorm
!
      Real
     &           rnorm_buf
!
!     %-----------------------%
!     | Local Array Arguments |
!     %-----------------------%
!
      Real
     &           xtemp(2)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   saxpy, scopy, sscal, sgemv, psgetv0, psvout, psmout,
     &           slascl, pivout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real
     &           sdot, psnorm2, pslamch10
      external   sdot, psnorm2, pslamch10
!
!     %-----------------%
!     | Data statements |
!     %-----------------%
!
!
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (aitr_first) then
         aitr_first = .false.
!
!        %--------------------------------%
!        | safmin = safe minimum is such  |
!        | that 1/sfmin does not overflow |
!        %--------------------------------%
!
         safmin = pslamch10(comm,'safmin')
      end if
!
      if (ido .eq. 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call arscnd (t0)
         msglvl = msaitr
!
!        %------------------------------%
!        | Initial call to this routine |
!        %------------------------------%
!
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
!
!        %--------------------------------%
!        | Pointer to the current step of |
!        | the factorization to build     |
!        %--------------------------------%
!
         j      = k + 1
!
!        %------------------------------------------%
!        | Pointers used for reverse communication  |
!        | when using WORKD.                        |
!        %------------------------------------------%
!
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
!
!     %-------------------------------------------------%
!     | When in reverse communication mode one of:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | will be .true.                                  |
!     | STEP3: return from computing OP*v_{j}.          |
!     | STEP4: return from computing B-norm of OP*v_{j} |
!     | ORTH1: return from computing B-norm of r_{j+1}  |
!     | ORTH2: return from computing B-norm of          |
!     |        correction to the residual vector.       |
!     | RSTART: return from OP computations needed by   |
!     |         psgetv0.                                |
!     %-------------------------------------------------%
!
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
!
!     %------------------------------%
!     | Else this is the first step. |
!     %------------------------------%
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%
!
 1000 continue
!
         if (msglvl .gt. 2) then
            call pivout (comm, logfil, 1, [j], ndigit,
     &                  '_saitr: generating Arnoldi vector no.')
            call psvout (comm, logfil, 1, [rnorm], ndigit,
     &                  '_saitr: B-norm of the current residual =')
         end if
!
!        %---------------------------------------------------------%
!        | Check for exact zero. Equivalent to determining whether a |
!        | j-step Arnoldi factorization is present.                |
!        %---------------------------------------------------------%
!
         if (rnorm .gt. zero) go to 40
!
!           %---------------------------------------------------%
!           | Invariant subspace found, generate a new starting |
!           | vector which is orthogonal to the current Arnoldi |
!           | basis and continue the iteration.                 |
!           %---------------------------------------------------%
!
            if (msglvl .gt. 0) then
               call pivout (comm, logfil, 1, [j], ndigit,
     &                     '_saitr: ****** restart at step ******')
            end if
!
!           %---------------------------------------------%
!           | ITRY is the loop variable that controls the |
!           | maximum amount of times that a restart is   |
!           | attempted. NRSTRT is used by stat.h         |
!           %---------------------------------------------%
!
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
!
!           %--------------------------------------%
!           | If in reverse communication mode and |
!           | RSTART = .true. flow returns here.   |
!           %--------------------------------------%
!
            call psgetv0 (comm, ido, bmat, itry, .false., n, j, v, ldv,
     &                   resid, rnorm, ipntr, workd, workl, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
!
!              %------------------------------------------------%
!              | Give up after several restart attempts.        |
!              | Set INFO to the size of the invariant subspace |
!              | which spans OP and exit.                       |
!              %------------------------------------------------%
!
               info = j - 1
               call arscnd (t1)
               tsaitr = tsaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
!
   40    continue
!
!        %---------------------------------------------------------%
!        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
!        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
!        | when reciprocating a small RNORM, test against lower    |
!        | machine bound.                                          |
!        %---------------------------------------------------------%
!
         call scopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. safmin) then
             temp1 = one / rnorm
             call sscal (n, temp1, v(1,j), 1)
             call sscal (n, temp1, workd(ipj), 1)
         else
!
!            %-----------------------------------------%
!            | To scale both v_{j} and p_{j} carefully |
!            | use LAPACK routine SLASCL               |
!            %-----------------------------------------%
!
             call slascl ('General', i, i, rnorm, one, n, 1,
     &                    v(1,j), n, infol)
             call slascl ('General', i, i, rnorm, one, n, 1,
     &                    workd(ipj), n, infol)
         end if
!
!        %------------------------------------------------------%
!        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
!        | Note that this is not quite yet r_{j}. See STEP 4    |
!        %------------------------------------------------------%
!
         step3 = .true.
         nopx  = nopx + 1
         call arscnd (t2)
         call scopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
!
!        %-----------------------------------%
!        | Exit in order to compute OP*v_{j} |
!        %-----------------------------------%
!
         go to 9000
   50    continue
!
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
!        %-----------------------------------%
!
         call arscnd (t3)
         tmvopx = tmvopx + (t3 - t2)
!
         step3 = .false.
!
!        %------------------------------------------%
!        | Put another copy of OP*v_{j} into RESID. |
!        %------------------------------------------%
!
         call scopy (n, workd(irj), 1, resid, 1)
!
!        %-------------------------------------------%
!        | STEP 4:  Finish extending the symmetric   |
!        |          Arnoldi to length j. If MODE = 2 |
!        |          then B*OP = B*inv(B)*A = A and   |
!        |          we don't need to compute B*OP.   |
!        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
!        | assumed to have A*v_{j}.                  |
!        %-------------------------------------------%
!
         if (mode .eq. 2) go to 65
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %-------------------------------------%
!           | Exit in order to compute B*OP*v_{j} |
!           %-------------------------------------%
!
            go to 9000
         else if (bmat .eq. 'I') then
              call scopy(n, resid, 1 , workd(ipj), 1)
         end if
   60    continue
!
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
!        %-----------------------------------%
!
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         step4 = .false.
!
!        %-------------------------------------%
!        | The following is needed for STEP 5. |
!        | Compute the B-norm of OP*v_{j}.     |
!        %-------------------------------------%
!
   65    continue
         if (mode .eq. 2) then
!
!           %----------------------------------%
!           | Note that the B-norm of OP*v_{j} |
!           | is the inv(B)-norm of A*v_{j}.   |
!           %----------------------------------%
!
            rnorm_buf = sdot (n, resid, 1, workd(ivj), 1)
            call MPI_ALLREDUCE( [rnorm_buf], wnorm, 1,
     &           MPI_REAL, MPI_SUM, comm, ierr )
            wnorm(1) = sqrt(abs(wnorm(1)))
         else if (bmat .eq. 'G') then
            rnorm_buf = sdot (n, resid, 1, workd(ipj), 1)
            call MPI_ALLREDUCE( [rnorm_buf], wnorm, 1,
     &           MPI_REAL, MPI_SUM, comm, ierr )
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm(1) = psnorm2( comm, n, resid, 1 )
         end if
!
!        %-----------------------------------------%
!        | Compute the j-th residual corresponding |
!        | to the j step factorization.            |
!        | Use Classical Gram Schmidt and compute: |
!        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
!        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
!        %-----------------------------------------%
!
!
!        %------------------------------------------%
!        | Compute the j Fourier coefficients w_{j} |
!        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
!        %------------------------------------------%
!
         if (mode .ne. 2 ) then
            call sgemv('T', n, j, one, v, ldv, workd(ipj), 1, zero,
     &                  workl(j+1), 1)
            call MPI_ALLREDUCE( workl(j+1), workl(1), j,
     &                  MPI_REAL, MPI_SUM, comm, ierr)
         else if (mode .eq. 2) then
            call sgemv('T', n, j, one, v, ldv, workd(ivj), 1, zero,
     &                  workl(j+1), 1)
            call MPI_ALLREDUCE( workl(j+1), workl(1), j,
     &                  MPI_REAL, MPI_SUM, comm, ierr)
         end if
!
!        %--------------------------------------%
!        | Orthgonalize r_{j} against V_{j}.    |
!        | RESID contains OP*v_{j}. See STEP 3. |
!        %--------------------------------------%
!
         call sgemv('N', n, j, -one, v, ldv, workl(1), 1, one,
     &               resid, 1)
!
!        %--------------------------------------%
!        | Extend H to have j rows and columns. |
!        %--------------------------------------%
!
         h(j,2) = workl(j)
         if (j .eq. 1  .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         end if
         call arscnd (t4)
!
         orth1 = .true.
         iter  = 0
!
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*r_{j} |
!           %----------------------------------%
!
            go to 9000
         else if (bmat .eq. 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH1 = .true. |
!        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
!        %---------------------------------------------------%
!
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         orth1 = .false.
!
!        %------------------------------%
!        | Compute the B-norm of r_{j}. |
!        %------------------------------%
!
         if (bmat .eq. 'G') then
            rnorm_buf = sdot (n, resid, 1, workd(ipj), 1)
            call MPI_ALLREDUCE( [rnorm_buf], temp2, 1,
     &           MPI_REAL, MPI_SUM, comm, ierr )
            rnorm = sqrt(abs(temp2(1)))
         else if (bmat .eq. 'I') then
            rnorm = psnorm2( comm, n, resid, 1 )
         end if
!
!        %-----------------------------------------------------------%
!        | STEP 5: Re-orthogonalization / Iterative refinement phase |
!        | Maximum NITER_ITREF tries.                                |
!        |                                                           |
!        |          s      = V_{j}^T * B * r_{j}                     |
!        |          r_{j}  = r_{j} - V_{j}*s                         |
!        |          alphaj = alphaj + s_{j}                          |
!        |                                                           |
!        | The stopping criteria used for iterative refinement is    |
!        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
!        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
!        | Determine if we need to correct the residual. The goal is |
!        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
!        %-----------------------------------------------------------%
!
         if (rnorm .gt. 0.717*wnorm(1)) go to 100
         nrorth = nrorth + 1
!
!        %---------------------------------------------------%
!        | Enter the Iterative refinement phase. If further  |
!        | refinement is necessary, loop back here. The loop |
!        | variable is ITER. Perform a step of Classical     |
!        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
!        %---------------------------------------------------%
!
   80    continue
!
         if (msglvl .gt. 2) then
            xtemp(1) = wnorm(1)
            xtemp(2) = rnorm
            call psvout (comm, logfil, 2, xtemp, ndigit,
     &           '_naitr: re-orthonalization ; wnorm and rnorm are')
         end if
!
!        %----------------------------------------------------%
!        | Compute V_{j}^T * B * r_{j}.                       |
!        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
         call sgemv ('T', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, workl(j+1), 1)
         call MPI_ALLREDUCE( workl(j+1), workl(1), j,
     &               MPI_REAL, MPI_SUM, comm, ierr)
!
!        %----------------------------------------------%
!        | Compute the correction to the residual:      |
!        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
!        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
!        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
!        | H(j,j) is updated.                           |
!        %----------------------------------------------%
!
         call sgemv ('N', n, j, -one, v, ldv, workl(1), 1,
     &               one, resid, 1)
!
         if (j .eq. 1  .or.  rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workl(j)
!
         orth2 = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %-----------------------------------%
!           | Exit in order to compute B*r_{j}. |
!           | r_{j} is the corrected residual.  |
!           %-----------------------------------%
!
            go to 9000
         else if (bmat .eq. 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH2 = .true. |
!        %---------------------------------------------------%
!
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
!        %-----------------------------------------------------%
!        | Compute the B-norm of the corrected residual r_{j}. |
!        %-----------------------------------------------------%
!
         if (bmat .eq. 'G') then
           rnorm_buf = sdot (n, resid, 1, workd(ipj), 1)
           call MPI_ALLREDUCE( [rnorm_buf], temp2, 1,
     &          MPI_REAL, MPI_SUM, comm, ierr )
           rnorm1 = sqrt(abs(temp2(1)))
         else if (bmat .eq. 'I') then
           rnorm1 = psnorm2( comm, n, resid, 1 )
         end if
!
         if (msglvl .gt. 0 .and. iter .gt. 0) then
            call pivout (comm, logfil, 1, [j], ndigit,
     &           '_naitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call psvout (comm, logfil, 2, xtemp, ndigit,
     &           '_naitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
!
!        %-----------------------------------------%
!        | Determine if we need to perform another |
!        | step of re-orthogonalization.           |
!        %-----------------------------------------%
!
         if (rnorm1 .gt. 0.717*rnorm) then
!
!           %--------------------------------%
!           | No need for further refinement |
!           %--------------------------------%
!
            rnorm = rnorm1
!
         else
!
!           %-------------------------------------------%
!           | Another step of iterative refinement step |
!           | is required. NITREF is used by stat.h     |
!           %-------------------------------------------%
!
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
!
!           %-------------------------------------------------%
!           | Otherwise RESID is numerically in the span of V |
!           %-------------------------------------------------%
!
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
!
!        %----------------------------------------------%
!        | Branch here directly if iterative refinement |
!        | wasn't necessary or after at most NITER_REF  |
!        | steps of iterative refinement.               |
!        %----------------------------------------------%
!
  100    continue
!
         rstart = .false.
         orth2  = .false.
!
         call arscnd (t5)
         titref = titref + (t5 - t4)
!
!        %----------------------------------------------------------%
!        | Make sure the last off-diagonal element is non negative  |
!        | If not perform a similarity transformation on H(1:j,1:j) |
!        | and scale v(:,j) by -1.                                  |
!        %----------------------------------------------------------%
!
         if (h(j,1) .lt. zero) then
            h(j,1) = -h(j,1)
            if ( j .lt. k+np) then
               call sscal(n, -one, v(1,j+1), 1)
            else
               call sscal(n, -one, resid, 1)
            end if
         end if
!
!        %------------------------------------%
!        | STEP 6: Update  j = j+1;  Continue |
!        %------------------------------------%
!
         j = j + 1
         if (j .gt. k+np) then
            call arscnd (t1)
            tsaitr = tsaitr + (t1 - t0)
            ido = 99
!
            if (msglvl .gt. 1) then
               call psvout (comm, logfil, k+np, h(1,2), ndigit,
     &         '_saitr: main diagonal of matrix H of step K+NP.')
               if (k+np .gt. 1) then
               call psvout (comm, logfil, k+np-1, h(2,1), ndigit,
     &         '_saitr: sub diagonal of matrix H of step K+NP.')
               end if
            end if
!
            go to 9000
         end if
!
!        %--------------------------------------------------------%
!        | Loop back to extend the factorization by another step. |
!        %--------------------------------------------------------%
!
      go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 9000 continue
      return
!
!     %----------------%
!     | End of pssaitr |
!     %----------------%
!
      end
