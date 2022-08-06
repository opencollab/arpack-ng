!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: psgetv0
!
! Message Passing Layer: MPI
!
!\Description:
!  Generate a random initial residual vector for the Arnoldi process.
!  Force the residual vector to be in the range of the operator OP.
!
!\Usage:
!  call psgetv0
!     ( COMM, IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
!       IPNTR, WORKD, WORKL, IERR )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to psgetv0.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B in the (generalized)
!          eigenvalue problem A*x = lambda*B*x.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  ITRY    Integer.  (INPUT)
!          ITRY counts the number of times that psgetv0 is called.
!          It should be set to 1 on the initial call to psgetv0.
!
!  INITV   Logical variable.  (INPUT)
!          .TRUE.  => the initial residual vector is given in RESID.
!          .FALSE. => generate a random initial residual vector.
!
!  N       Integer.  (INPUT)
!          Dimension of the problem.
!
!  J       Integer.  (INPUT)
!          Index of the residual vector to be generated, with respect to
!          the Arnoldi process.  J > 1 in case of a "restart".
!
!  V       Real N by J array.  (INPUT)
!          The first J-1 columns of V contain the current Arnoldi basis
!          if this is a "restart".
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          Initial residual vector to be generated.  If RESID is
!          provided, force RESID into the range of the operator OP.
!
!  RNORM   Real scalar.  (OUTPUT)
!          B-norm of the generated residual.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!
!  WORKD   Real work array of length 2*N.  (REVERSE COMMUNICATION).
!          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
!
!  WORKL   Real work space used for Gram Schmidt orthogonalization
!
!  IERR    Integer.  (OUTPUT)
!          =  0: Normal exit.
!          = -1: Cannot generate a nontrivial restarted residual vector
!                in the range of the operator OP.
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
!     arscnd   ARPACK utility routine for timing.
!     psvout   Parallel ARPACK utility routine for vector output.
!     pslarnv  Parallel wrapper for LAPACK routine slarnv (generates a random vector).
!     sgemv    Level 2 BLAS routine for matrix vector multiplication.
!     scopy    Level 1 BLAS that copies one vector to another.
!     sdot     Level 1 BLAS that computes the scalar product of two vectors.
!     psnorm2  Parallel version of  Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              Cray Research, Inc. &
!     Dept. of Computational &     CRPC / Rice University
!     Applied Mathematics          Houston, Texas
!     Rice University
!     Houston, Texas
!
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Serial Code FILE: getv0.F   SID: 2.3
!
!\SCCS Information:
! FILE: getv0.F   SID: 1.4   DATE OF SID: 3/19/97
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine psgetv0&
         ( comm, ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm,&
           ipntr, workd, workl, ierr )
!
#ifdef HAVE_MPI_ICB
      use :: mpi_f08
#else
#include "mpif.h"
#endif

!
!     %---------------%
!     | MPI Variables |
!     %---------------%
!
#ifdef HAVE_MPI_ICB
      type(MPI_Comm) comm
#else
      integer    comm
#endif
      integer*4  cnt, ierr

!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
#include "debugF90.h"
#include "statF90.h"
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1
      logical    initv
      integer    ido, itry, j, ldv, n
      Real&
                 rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Real&
                 resid(n), v(ldv,j), workd(2*n), workl(2*j), buf2(1)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real&
                 one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    first, orth
      integer    idist, iseed(4), iter, msglvl, jj
      Real&
                 rnorm0
      save       first, iseed, iter, msglvl, orth, rnorm0
!
      Real&
                 rnorm_buf
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   pslarnv, psvout, scopy, sgemv, arscnd
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real&
                 sdot, psnorm2
      external   sdot, psnorm2
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
!
!     %-----------------------------------%
!     | Initialize the seed of the LAPACK |
!     | random number generator           |
!     %-----------------------------------%
!
      iseed(1) = 1
      iseed(2) = 3
      iseed(3) = 5
      iseed(4) = 7
!
      if (ido .eq.  0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call arscnd (t0)
         msglvl = mgetv0
!
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
!
!        %-----------------------------------------------------%
!        | Possibly generate a random starting vector in RESID |
!        | Use a LAPACK random number generator used by the    |
!        | matrix generation routines.                         |
!        |    idist = 1: uniform (0,1)  distribution;          |
!        |    idist = 2: uniform (-1,1) distribution;          |
!        |    idist = 3: normal  (0,1)  distribution;          |
!        %-----------------------------------------------------%
!
         if (.not.initv) then
            idist = 2
            call pslarnv (comm, idist, iseed, n, resid)
         end if
!
!        %----------------------------------------------------------%
!        | Force the starting vector into the range of OP to handle |
!        | the generalized problem when B is possibly (singular).   |
!        %----------------------------------------------------------%
!
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call scopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
!
!     %-----------------------------------------%
!     | Back from computing OP*(initial-vector) |
!     %-----------------------------------------%
!
      if (first) go to 20
!
!     %-----------------------------------------------%
!     | Back from computing B*(orthogonalized-vector) |
!     %-----------------------------------------------%
!
      if (orth)  go to 40
!
      call arscnd (t3)
      tmvopx = tmvopx + (t3 - t2)
!
!     %------------------------------------------------------%
!     | Starting vector is now in the range of OP; r = OP*r; |
!     | Compute B-norm of starting vector.                   |
!     %------------------------------------------------------%
!
      call arscnd (t2)
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call scopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call scopy (n, resid, 1, workd, 1)
      end if
!
   20 continue
!
      if (bmat .eq. 'G') then
         call arscnd (t3)
         tmvbx = tmvbx + (t3 - t2)
      endif
!
      first = .FALSE.
      if (bmat .eq. 'G') then
          rnorm_buf = sdot (n, resid, 1, workd, 1)
          cnt = 1
          call MPI_ALLREDUCE( [rnorm_buf], buf2, cnt,&
                MPI_REAL, MPI_SUM, comm, ierr )
          rnorm0 = sqrt(abs(buf2(1)))
      else if (bmat .eq. 'I') then
          rnorm0 = psnorm2( comm, n, resid, 1 )
      end if
      rnorm  = rnorm0
!
!     %---------------------------------------------%
!     | Exit if this is the very first Arnoldi step |
!     %---------------------------------------------%
!
      if (j .eq. 1) go to 50
!
!     %----------------------------------------------------------------
!     | Otherwise need to B-orthogonalize the starting vector against |
!     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
!     | This is the case where an invariant subspace is encountered   |
!     | in the middle of the Arnoldi factorization.                   |
!     |                                                               |
!     |       s = V^{T}*B*r;   r = r - V*s;                           |
!     |                                                               |
!     | Stopping criteria used for iter. ref. is discussed in         |
!     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
!     %---------------------------------------------------------------%
!
      orth = .TRUE.
   30 continue
!
      call sgemv ('T', n, j-1, one, v, ldv, workd, 1,&
                  zero, workl(j+1), 1)
      cnt = j-1
      call MPI_ALLREDUCE( workl(j+1), workl, cnt,&
                          MPI_REAL, MPI_SUM, comm, ierr)
      call sgemv ('N', n, j-1, -one, v, ldv, workl, 1,&
                  one, resid, 1)
!
!     %----------------------------------------------------------%
!     | Compute the B-norm of the orthogonalized starting vector |
!     %----------------------------------------------------------%
!
      call arscnd (t2)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call scopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call scopy (n, resid, 1, workd, 1)
      end if
!
   40 continue
!
      if (bmat .eq. 'G') then
         call arscnd (t3)
         tmvbx = tmvbx + (t3 - t2)
      endif
!
      if (bmat .eq. 'G') then
         rnorm_buf = sdot (n, resid, 1, workd, 1)
         cnt = 1
         call MPI_ALLREDUCE( [rnorm_buf], buf2, cnt,&
                  MPI_REAL, MPI_SUM, comm, ierr )
         rnorm = sqrt(abs(buf2(1)))
      else if (bmat .eq. 'I') then
         rnorm = psnorm2( comm, n, resid, 1 )
      end if
!
!     %--------------------------------------%
!     | Check for further orthogonalization. |
!     %--------------------------------------%
!
      if (msglvl .gt. 2) then
          call psvout (comm, logfil, 1, [rnorm0], ndigit,&
                      '_getv0: re-orthonalization ; rnorm0 is')
          call psvout (comm, logfil, 1, [rnorm], ndigit,&
                      '_getv0: re-orthonalization ; rnorm is')
      end if
!
      if (rnorm .gt. 0.717*rnorm0) go to 50
!
      iter = iter + 1
      if (iter .le. 1) then
!
!        %-----------------------------------%
!        | Perform iterative refinement step |
!        %-----------------------------------%
!
         rnorm0 = rnorm
         go to 30
      else
!
!        %------------------------------------%
!        | Iterative refinement step "failed" |
!        %------------------------------------%
!
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = zero
         ierr = -1
      end if
!
   50 continue
!
      if (msglvl .gt. 0) then
         call psvout (comm, logfil, 1, [rnorm], ndigit,&
              '_getv0: B-norm of initial / restarted starting vector')
      end if
      if (msglvl .gt. 2) then
         call psvout (comm, logfil, n, resid, ndigit,&
              '_getv0: initial / restarted starting vector')
      end if
      ido = 99
!
      call arscnd (t1)
      tgetv0 = tgetv0 + (t1 - t0)
!
 9000 continue
      return
!
!     %----------------%
!     | End of psgetv0 |
!     %----------------%
!
      end
