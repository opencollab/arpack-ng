!\BeginDoc
!
!\Name: pznaitr
!
! Message Passing Layer: MPI
!
!\Description:
!  Reverse communication interface for applying NP additional steps to
!  a K step nonsymmetric Arnoldi factorization.
!
!  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
!
!          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
!
!  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
!
!          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
!
!  where OP and B are as in pznaupd.  The B-norm of r_{k+p} is also
!  computed and returned.
!
!\Usage:
!  call pznaitr
!     ( COMM, IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH,
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
!          vector B * Q is already available and do not need to be
!          recomputed in forming OP * Q.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.  See pznaupd.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  K       Integer.  (INPUT)
!          Current size of V and H.
!
!  NP      Integer.  (INPUT)
!          Number of additional Arnoldi steps to take.
!
!  NB      Integer.  (INPUT)
!          Blocksize to be used in the recurrence.
!          Only work for NB = 1 right now.  The goal is to have a
!          program that implement both the block and non-block method.
!
!  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
!          On INPUT:  RESID contains the residual vector r_{k}.
!          On OUTPUT: RESID contains the residual vector r_{k+p}.
!
!  RNORM   Double precision scalar.  (INPUT/OUTPUT)
!          B-norm of the starting residual on input.
!          B-norm of the updated residual r_{k+p} on output.
!
!  V       Complex*16 N by K+NP array.  (INPUT/OUTPUT)
!          On INPUT:  V contains the Arnoldi vectors in the first K
!          columns.
!          On OUTPUT: V contains the new NP Arnoldi vectors in the next
!          NP columns.  The first K columns are unchanged.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Complex*16 (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
!          H is used to store the generated upper Hessenberg matrix.
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
!  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The calling program should not
!          use WORKD as temporary workspace during the iteration !!!!!!
!          On input, WORKD(1:N) = B*RESID and is used to save some
!          computation at the first step.
!
!  WORKL   Complex*16 work space used for Gram Schmidt orthogonalization
!
!  INFO    Integer.  (OUTPUT)
!          = 0: Normal exit.
!          > 0: Size of the spanning invariant subspace of OP found.
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
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     pzgetv0   Parallel ARPACK routine to generate the initial vector.
!     pivout    Parallel ARPACK utility routine that prints integers.
!     arscnd    ARPACK utility routine for timing.
!     pzmout    Parallel ARPACK utility routine that prints matrices
!     pzvout    Parallel ARPACK utility routine that prints vectors.
!     zlanhs    LAPACK routine that computes various norms of a matrix.
!     zlascl    LAPACK routine for careful scaling of a matrix.
!     dlabad    LAPACK routine for defining the underflow and overflow
!               limits
!     pdlamch10   ScaLAPACK routine that determines machine constants.
!     dlapy2    LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     zgemv     Level 2 BLAS routine for matrix vector multiplication.
!     zaxpy     Level 1 BLAS that computes a vector triad.
!     zcopy     Level 1 BLAS that copies one vector to another .
!     zdotc     Level 1 BLAS that computes the scalar product of
!               two vectors.
!     zscal     Level 1 BLAS that scales a vector.
!     zdscal    Level 1 BLAS that scales a complex vector by a real number.
!     pdznorm2  Parallel version of Level 1 BLAS that computes the
!               norm of a vector.
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
!     Starting Point: Complex Code FILE: naitr.F   SID: 2.1
!
!\SCCS Information:
! FILE: naitr.F   SID: 1.3   DATE OF SID: 3/19/97
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
!        ( At present tol is zero )
!        if ( restart ) generate a new starting vector.
!     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
!        p_{j} = p_{j}/betaj
!     3) r_{j} = OP*v_{j} where OP is defined as in pznaupd
!        For shift-invert mode p_{j} = B*v_{j} is already available.
!        wnorm = || OP*v_{j} ||
!     4) Compute the j-th step residual vector.
!        w_{j} =  V_{j}^T * B * OP * v_{j}
!        r_{j} =  OP*v_{j} - V_{j} * w_{j}
!        H(:,j) = w_{j};
!        H(j,j-1) = rnorm
!        rnorm = || r_(j) ||
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
      subroutine pznaitr
     &   (comm, ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh,
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
      integer    ido, info, k, ldh, ldv, n, nb, np
      Double precision
     &           rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Complex*16
     &           h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n),
     &           workl(2*ldh)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex*16
     &           one, zero
      Double precision
     &           rone, rzero
      parameter (one = (1.0, 0.0), zero = (0.0, 0.0),
     &           rone = 1.0, rzero = 0.0)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      Double precision
     &           rtemp(2)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      logical    orth1, orth2, rstart, step3, step4
      integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j, msglvl,
     &           jj
      Double precision
     &           ovfl, smlnum, tst1, ulp, unfl, betaj,
     &           temp1, rnorm1, wnorm
      Complex*16
     &           cnorm
!
      save       orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl, ovfl,
     &           betaj, rnorm1, smlnum, ulp, unfl, wnorm
!
      Complex*16
     &           cnorm_buf, buf2(1)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   zaxpy, zcopy, zscal, zgemv, pzgetv0, dlabad,
     &           zdscal, pzvout, pzmout, pivout, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Complex*16
     &           zzdotc
      Double precision
     &           pdlamch10, pdznorm2, zlanhs, dlapy2
      external   zzdotc, pdznorm2, zlanhs, pdlamch10, dlapy2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic  dimag, dble, max, sqrt
!
!     %-----------------%
!     | Data statements |
!     %-----------------%
!
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (aitr_first) then
!
!        %-----------------------------------------%
!        | Set machine-dependent constants for the |
!        | the splitting and deflation criterion.  |
!        | If norm(H) <= sqrt(OVFL),               |
!        | overflow should not occur.              |
!        | REFERENCE: LAPACK subroutine zlahqr     |
!        %-----------------------------------------%
!
         unfl = pdlamch10(comm,  'safe minimum' )
         ovfl = dble(one / unfl)
         call dlabad( unfl, ovfl )
         ulp = pdlamch10( comm, 'precision' )
         smlnum = unfl*( n / ulp )
         aitr_first = .false.
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
         msglvl = mcaitr
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
         j      = k + 1
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
!
!     %-------------------------------------------------%
!     | When in reverse communication mode one of:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | will be .true. when ....                        |
!     | STEP3: return from computing OP*v_{j}.          |
!     | STEP4: return from computing B-norm of OP*v_{j} |
!     | ORTH1: return from computing B-norm of r_{j+1}  |
!     | ORTH2: return from computing B-norm of          |
!     |        correction to the residual vector.       |
!     | RSTART: return from OP computations needed by   |
!     |         pzgetv0.                                |
!     %-------------------------------------------------%
!
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
!
!     %-----------------------------%
!     | Else this is the first step |
!     %-----------------------------%
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%

 1000 continue
!
         if (msglvl .gt. 1) then
            call pivout (comm, logfil, 1, [j], ndigit,
     &                  '_naitr: generating Arnoldi vector number')
            call pdvout (comm, logfil, 1, [rnorm], ndigit,
     &                  '_naitr: B-norm of the current residual is')
         end if
!
!        %---------------------------------------------------%
!        | STEP 1: Check if the B norm of j-th residual      |
!        | vector is zero. Equivalent to determine whether   |
!        | an exact j-step Arnoldi factorization is present. |
!        %---------------------------------------------------%
!
         betaj = rnorm
         if (rnorm .gt. rzero) go to 40
!
!           %---------------------------------------------------%
!           | Invariant subspace found, generate a new starting |
!           | vector which is orthogonal to the current Arnoldi |
!           | basis and continue the iteration.                 |
!           %---------------------------------------------------%
!
            if (msglvl .gt. 0) then
               call pivout (comm, logfil, 1, [j], ndigit,
     &                     '_naitr: ****** RESTART AT STEP ******')
            end if
!
!           %---------------------------------------------%
!           | ITRY is the loop variable that controls the |
!           | maximum amount of times that a restart is   |
!           | attempted. NRSTRT is used by stat.h         |
!           %---------------------------------------------%
!
            betaj  = rzero
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
            call pzgetv0 (comm, ido, bmat, itry, .false., n, j, v, ldv,
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
               tcaitr = tcaitr + (t1 - t0)
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
         call zcopy (n, resid, 1, v(1,j), 1)
         if ( rnorm .ge. unfl) then
             temp1 = rone / rnorm
             call zdscal (n, temp1, v(1,j), 1)
             call zdscal (n, temp1, workd(ipj), 1)
         else
!
!            %-----------------------------------------%
!            | To scale both v_{j} and p_{j} carefully |
!            | use LAPACK routine zlascl               |
!            %-----------------------------------------%
!
             call zlascl ('General', i, i, rnorm, rone,
     &                    n, 1, v(1,j), n, infol)
             call zlascl ('General', i, i, rnorm, rone,
     &                    n, 1, workd(ipj), n, infol)
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
         call zcopy (n, v(1,j), 1, workd(ivj), 1)
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
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
!        | if step3 = .true.                |
!        %----------------------------------%
!
         call arscnd (t3)
         tmvopx = tmvopx + (t3 - t2)

         step3 = .false.
!
!        %------------------------------------------%
!        | Put another copy of OP*v_{j} into RESID. |
!        %------------------------------------------%
!
         call zcopy (n, workd(irj), 1, resid, 1)
!
!        %---------------------------------------%
!        | STEP 4:  Finish extending the Arnoldi |
!        |          factorization to length j.   |
!        %---------------------------------------%
!
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
            call zcopy (n, resid, 1, workd(ipj), 1)
         end if
   60    continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
!        | if step4 = .true.                |
!        %----------------------------------%
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
         if (bmat .eq. 'G') then
             cnorm_buf = zzdotc (n, resid, 1, workd(ipj), 1)
            call MPI_ALLREDUCE( [cnorm_buf], buf2, 1,
     &           MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr )
            cnorm = buf2(1)
            wnorm = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
             wnorm = pdznorm2(comm, n, resid, 1)
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
         call zgemv ('C', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, workl, 1)
         call MPI_ALLREDUCE( workl, h(1,j), j,
     &               MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
!
!        %--------------------------------------%
!        | Orthogonalize r_{j} against V_{j}.   |
!        | RESID contains OP*v_{j}. See STEP 3. |
!        %--------------------------------------%
!
         call zgemv ('N', n, j, -one, v, ldv, h(1,j), 1,
     &               one, resid, 1)
!
         if (j .gt. 1) h(j,j-1) = dcmplx(betaj, rzero)
!
         call arscnd (t4)
!
         orth1 = .true.
!
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy (n, resid, 1, workd(irj), 1)
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
            call zcopy (n, resid, 1, workd(ipj), 1)
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
            cnorm_buf = zzdotc (n, resid, 1, workd(ipj), 1)
            call MPI_ALLREDUCE( [cnorm_buf], buf2, 1,
     &           MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr )
            cnorm = buf2(1)
            rnorm = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
            rnorm = pdznorm2(comm, n, resid, 1)
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
!        | The following test determines whether the sine of the     |
!        | angle between  OP*x and the computed residual is less     |
!        | than or equal to 0.717.                                   |
!        %-----------------------------------------------------------%
!
         if ( rnorm .gt. 0.717*wnorm ) go to 100
!
         iter  = 0
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
            rtemp(1) = wnorm
            rtemp(2) = rnorm
            call pdvout (comm, logfil, 2, rtemp, ndigit,
     &      '_naitr: re-orthogonalization; wnorm and rnorm are')
            call pzvout (comm, logfil, j, h(1,j), ndigit,
     &                  '_naitr: j-th column of H')
         end if
!
!        %----------------------------------------------------%
!        | Compute V_{j}^T * B * r_{j}.                       |
!        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
         call zgemv ('C', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, workl(j+1), 1)
         call MPI_ALLREDUCE( workl(j+1), workl(1), j,
     &               MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
!
!        %---------------------------------------------%
!        | Compute the correction to the residual:     |
!        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
!        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
!        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
!        %---------------------------------------------%
!
         call zgemv ('N', n, j, -one, v, ldv, workl(1), 1,
     &               one, resid, 1)
         call zaxpy (j, one, workl(1), 1, h(1,j), 1)
!
         orth2 = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call zcopy (n, resid, 1, workd(irj), 1)
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
            call zcopy (n, resid, 1, workd(ipj), 1)
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
             cnorm_buf = zzdotc (n, resid, 1, workd(ipj), 1)
            call MPI_ALLREDUCE( [cnorm_buf], buf2, 1,
     &           MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr )
            cnorm = buf2(1)
            rnorm1 = sqrt( dlapy2(dble(cnorm),dimag(cnorm)) )
         else if (bmat .eq. 'I') then
             rnorm1 = pdznorm2(comm, n, resid, 1)
         end if
!
         if (msglvl .gt. 0 .and. iter .gt. 0 ) then
            call pivout (comm, logfil, 1, [j], ndigit,
     &           '_naitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                rtemp(1) = rnorm
                rtemp(2) = rnorm1
                call pdvout (comm, logfil, 2, rtemp, ndigit,
     &           '_naitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
!
!        %-----------------------------------------%
!        | Determine if we need to perform another |
!        | step of re-orthogonalization.           |
!        %-----------------------------------------%
!
         if ( rnorm1 .gt. 0.717*rnorm ) then
!
!           %---------------------------------------%
!           | No need for further refinement.       |
!           | The cosine of the angle between the   |
!           | corrected residual vector and the old |
!           | residual vector is greater than 0.717 |
!           | In other words the corrected residual |
!           | and the old residual vector share an  |
!           | angle of less than arcCOS(0.717)      |
!           %---------------------------------------%
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
            rnorm = rzero
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
!        %------------------------------------%
!        | STEP 6: Update  j = j+1;  Continue |
!        %------------------------------------%
!
         j = j + 1
         if (j .gt. k+np) then
            call arscnd (t1)
            tcaitr = tcaitr + (t1 - t0)
            ido = 99
            do 110 i = max(1,k), k+np-1
!
!              %--------------------------------------------%
!              | Check for splitting and deflation.         |
!              | Use a standard test as in the QR algorithm |
!              | REFERENCE: LAPACK subroutine zlahqr        |
!              %--------------------------------------------%
!
               tst1 = dlapy2(dble(h(i,i)),dimag(h(i,i)))
     &              + dlapy2(dble(h(i+1,i+1)), dimag(h(i+1,i+1)))
               if( tst1.eq.dble(zero) )
     &              tst1 = zlanhs( '1', k+np, h, ldh, workd(n+1) )
               if( dlapy2(dble(h(i+1,i)),dimag(h(i+1,i))) .le.
     &                    max( ulp*tst1, smlnum ) )
     &             h(i+1,i) = zero
 110        continue
!
            if (msglvl .gt. 2) then
               call pzmout (comm, logfil, k+np, k+np, h, ldh, ndigit,
     &          '_naitr: Final upper Hessenberg matrix H of order K+NP')
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
!     | End of pznaitr |
!     %----------------%
!
      end
