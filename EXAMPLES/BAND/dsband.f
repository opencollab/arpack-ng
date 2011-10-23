c \BeginDoc
c
c \Name: dsband
c
c \Description:
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  Matrices A and B are stored in LAPACK-style band form.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost 
c  of n*nev if both are requested (in this case a separate array Z must be 
c  supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are called Ritz values and Ritz vectors respectively.  They are referred 
c  to as such in the comments that follow.  The computed orthonormal basis 
c  for the invariant subspace corresponding to these Ritz values is referred 
c  to as a Lanczos basis.
c
c  dsband can be called with one of the following modes:
c
c  Mode 1:  A*x = lambda*x, A symmetric 
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 in DSAUPD)
c
c  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
c           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
c           ===> Shift-and-Invert mode
c
c  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite, 
c           KG symmetric indefinite
c           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
c           ===> Buckling mode
c
c  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
c           ===> Cayley transformed mode
c
c  The choice of mode must be specified in IPARAM(7) defined below.
c
c \Usage
c   call dsband
c      ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA, 
c        RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V, 
c        LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
c
c \Arguments
c
c  RVEC    Logical (INPUT)
c          Specifies whether Ritz vectors corresponding to the Ritz value 
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the associated Ritz vectors. 
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute all Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' , SELECT is not referenced.
c
c  D       Double precision array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd transformed to
c          those of the original eigensystem A*z = lambda*B*z. If 
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c          as the those of A*z = lambda*B*z.
c
c  Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Lanczos basis array V computed by DSAUPD.
c
c  LDZ     Integer.  (INPUT) 
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   Double precision  (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c 
c  N       Integer.  (INPUT) 
c          Dimension of the eigenproblem.  
c 
c  AB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set. 
c          The j-th column of M is stored in the j-th column of the
c          array AB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7) = 1
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, RFAC.
c
c  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
c          RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
c          It is not referenced when IPARAM(7) = 1.
c
c  KL      Integer. (INPUT)
c          Max(number of subdiagonals of A, number of subdiagonals of M)
c
c  KU      Integer. (OUTPUT)
c          Max(number of superdiagonals of A, number of superdiagonals of M)
c
c  WHICH   Character*2.  (INPUT)
c          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
c          the following.
c  
c            'LM' -> want the NEV eigenvalues of largest magnitude.
c            'SM' -> want the NEV eigenvalues of smallest magnitude.
c            'LA' -> want the NEV eigenvalues of largest REAL part.
c            'SA' -> want the NEV eigenvalues of smallest REAL part.
c            'BE' -> Compute NEV eigenvalues, half from each end of the 
c                    spectrum.  When NEV is odd, compute one more from 
c                    the high end than from the low end. 
c
c          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only. 
c          
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

c  NEV     Integer. (INPUT)
c          Number of eigenvalues of OP to be computed.
c   
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          Represents the dimension of the Lanczos basis constructed
c          by dsaupd for OP.
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Upon INPUT: the NCV columns of V contain the Lanczos basis 
c                      vectors as constructed by dsaupd for OP.
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
c                       represent the Ritz vectors that span the desired 
c                       invariant subspace.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Lanczos basis vector array V computed by dsaupd. In this case
c          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
c          the desired Ritz vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: 
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          It is set to 1 in this subroutine.  The user do not need
c          to set this parameter.
c          ------------------------------------------------------------
c          ISHIFT = 1: exact shifts with respect to the reduced 
c                      tridiagonal matrix T.  This is equivalent to 
c                      restarting the iteration with a starting vector 
c                      that is a linear combination of Ritz vectors 
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced. 
c
c          IPARAM(3) = MXITER
c          On INPUT:  max number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" eigenvalues.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used. 
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsband for the 
c          five modes available.
c
c          IPARAM(8) = NP
c          Not referenced.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
c
c WORKL    Double precision work array of length LWORKL.  (WORKSPACE)
c
c LWORKL   Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV.
c
c IWORK    Integer array of dimension at least N. (WORKSPACE)
c          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the 
c          factorization of M or (A-SIGMA*M).
c            
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 in DSAUPD.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informational error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -13: HOWMNY must be one of 'A' or 'P'
c          = -14: DSAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current
c                   Arnoldi factorization.
c
c \EndDoc
c
c------------------------------------------------------------------------
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
c     May 1995.
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dgbtrf  LAPACK band matrix factorization routine.
c     dgbtrs  LAPACK band linear system solve routine. 
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dgbmv   Level 2 BLAS that computes the band matrix vector product.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing 
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it.
c
c\Author    
c     Danny Sorensen
c     Richard Lehoucq
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: sband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
c
c\EndLib
c
c---------------------------------------------------------------------
c
      subroutine dsband( rvec, howmny, select, d, z, ldz, sigma, 
     &           n, ab, mb, lda, rfac, kl, ku, which, bmat, nev, 
     &           tol, resid, ncv, v, ldv, iparam, workd, workl, 
     &           lworkl, iwork, info)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c 
      character        which*2, bmat, howmny
      integer          n, lda, kl, ku, nev, ncv, ldv,
     &                 ldz, lworkl, info  
      Double precision
     &                 tol, sigma
      logical          rvec
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer          iparam(*), iwork(*)
      logical          select(*)
      Double precision
     &                 d(*), resid(*), v(ldv,*), z(ldz,*),
     &                 ab(lda,*), mb(lda,*), rfac(lda,*), 
     &                 workd(*), workl(*)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          ipntr(14)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer          ido, i, j, type, imid, itop, ibot, ierr
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &                  one, zero
      parameter        (one = 1.0, zero = 0.0)
c
c
c     %-----------------------------%
c     | LAPACK & BLAS routines used |
c     %-----------------------------%
c
      Double precision
     &                 ddot, dnrm2, dlapy2
      external         ddot, dcopy, dgbmv, dgbtrf, 
     &                 dgbtrs, dnrm2, dlapy2, dlacpy
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c     
c     %----------------------------------------------------------------%
c     | Set type of the problem to be solved. Check consistency        |
c     | between BMAT and IPARAM(7).                                    |
c     | type = 1 --> Solving standard problem in regular mode.         |
c     | type = 2 --> Solving standard problem in shift-invert mode.    | 
c     | type = 3 --> Solving generalized problem in regular mode.      |
c     | type = 4 --> Solving generalized problem in shift-invert mode. |
c     | type = 5 --> Solving generalized problem in Buckling mode.     |
c     | type = 6 --> Solving generalized problem in Cayley mode.       | 
c     %----------------------------------------------------------------%
c
      if ( iparam(7) .eq. 1 ) then
         type = 1
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'I') then
         type = 2
      else if ( iparam(7) .eq. 2 ) then
         type = 3
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
         type = 4 
      else if ( iparam(7) .eq. 4 ) then
         type = 5
      else if ( iparam(7) .eq. 5 ) then 
         type = 6
      else
         print*, ' '
         print*, 'BMAT is inconsistent with IPARAM(7).'
         print*, ' ' 
         go to 9000
      end if
c
c     %------------------------%
c     | Initialize the reverse |
c     | communication flag.    |
c     %------------------------%
c
      ido   = 0
c
c     %----------------%
c     | Exact shift is |
c     | used.          | 
c     %----------------%
c 
      iparam(1) = 1
c
c     %-----------------------------------%
c     | Both matrices A and M are stored  |
c     | between rows itop and ibot.  Imid |
c     | is the index of the row that      |
c     | stores the diagonal elements.     |
c     %-----------------------------------%
c
      itop = kl + 1
      imid = kl + ku + 1
      ibot = 2*kl + ku + 1
c
      if ( type .eq. 2 .or. type .eq. 6 .and. bmat .eq. 'I' ) then
c
c         %----------------------------------%
c         | Solving a standard eigenvalue    |
c         | problem in shift-invert or       |
c         | Cayley mode. Factor (A-sigma*I). |
c         %----------------------------------%
c
          call dlacpy ('A', ibot, n, ab, lda, rfac, lda )
          do 10 j = 1, n
             rfac(imid,j) =  ab(imid,j) - sigma
  10      continue
          call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr )
          if (ierr .ne. 0) then
             print*, ' ' 
             print*, ' _SBAND: Error with _gbtrf. '
             print*, ' '
             go to  9000
          end if
c
      else if ( type .eq. 3 ) then
c
c        %----------------------------------------------%
c        | Solving generalized eigenvalue problem in    |
c        | regular mode. Copy M to rfac and Call LAPACK |
c        | routine dgbtrf to factor M.                  |
c        %----------------------------------------------%
c
         call dlacpy ('A', ibot, n, mb, lda, rfac, lda )
         call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr) 
         if (ierr .ne. 0) then 
             print*, ' '
             print*,'_SBAND:  Error with _gbtrf.' 
             print*, ' ' 
             go to 9000 
         end if
c
      else if ( type .eq. 4 .or. type .eq. 5 .or. type .eq. 6 
     &         .and. bmat .eq. 'G' ) then
c
c        %-------------------------------------------%
c        | Solving generalized eigenvalue problem in |
c        | shift-invert, Buckling, or Cayley mode.   |
c        %-------------------------------------------%
c 
c        %-------------------------------------%
c        | Construct and factor (A - sigma*M). |
c        %-------------------------------------%
c
         do 60 j = 1,n
            do 50 i = itop, ibot 
               rfac(i,j) = ab(i,j) - sigma*mb(i,j)
  50        continue
  60     continue
c
         call dgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
         if ( ierr .ne. 0 )  then
             print*, ' '
             print*, '_SBAND: Error with _gbtrf.'
             print*, ' '
             go to 9000
         end if
c 
      end if 
c
c     %--------------------------------------------%
c     |  M A I N   L O O P (reverse communication) |
c     %--------------------------------------------%
c
  90  continue 
c
      call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
     &              info )
c
      if (ido .eq. -1) then
c
         if ( type .eq. 1) then
c
c           %----------------------------%
c           | Perform  y <--- OP*x = A*x |
c           %----------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( type .eq. 2 ) then
c
c           %----------------------------------%
c           |             Perform              |
c           | y <--- OP*x = inv[A-sigma*I]*x   |
c           | to force the starting vector     |
c           | into the range of OP.            |
c           %----------------------------------%
c
            call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                    iwork, workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' ' 
               print*, ' _SBAND: Error with _bgtrs. '
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 3 ) then
c
c           %-----------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x |
c           | to force the starting vector into | 
c           | the range of OP.                  |
c           %-----------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_SBAND: Error with sbgtrs.'
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 4 ) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x                      |
c           |           = inv[A-SIGMA*M]*M            | 
c           | to force the starting vector into the   |
c           | range of OP.                            |
c           %-----------------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' ' 
               print*, '_SBAND: Error with _gbtrs.'
               print*, ' ' 
               go to 9000
            end if
c
         else if ( type .eq. 5) then
c
c           %---------------------------------------% 
c           | Perform y <-- OP*x                    |
c           |    = inv[A-SIGMA*M]*A                 |
c           | to force the starting vector into the |
c           | range of OP.                          |
c           %---------------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
c
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _SBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 6 ) then
c
c           %---------------------------------------%
c           | Perform y <-- OP*x                    |
c           | = (inv[A-SIGMA*M])*(A+SIGMA*M)*x      | 
c           | to force the starting vector into the |
c           | range of OP.                          | 
c           %---------------------------------------%
c
            if ( bmat .eq. 'G' ) then
               call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                    ab(itop,1), lda, workd(ipntr(1)), 1, 
     &                    zero, workd(ipntr(2)), 1)
               call dgbmv('Notranspose', n, n, kl, ku, sigma, 
     &                    mb(itop,1), lda, workd(ipntr(1)), 1, 
     &                    one, workd(ipntr(2)), 1)
            else 
               call dcopy(n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                    lda, workd(ipntr(1)), 1, sigma, 
     &                    workd(ipntr(2)), 1)
            end if   
c
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
c
            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_SBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c
         end if
c
      else if (ido .eq. 1) then
c
         if ( type .eq. 1) then
c
c           %----------------------------%
c           | Perform  y <--- OP*x = A*x |
c           %----------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( type .eq. 2) then
c
c              %----------------------------------%
c              |             Perform              |
c              | y <--- OP*x = inv[A-sigma*I]*x.  |
c              %----------------------------------%
c
               call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                       iwork, workd(ipntr(2)), n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' '
                  print*, '_SBAND: Error with _gbtrs.' 
                  print*, ' '
                  go to 9000
               end if
c
         else if ( type .eq. 3 ) then
c
c           %-----------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x |
c           %-----------------------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1) 
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_SBAND: error with _bgtrs.'
               print*, ' ' 
               go to 9000
            end if
c
         else if ( type .eq. 4 ) then
c
c           %-------------------------------------%
c           | Perform y <-- inv(A-sigma*M)*(M*x). |
c           | (M*x) has been computed and stored  |
c           | in workd(ipntr(3)).                 |           
c           %-------------------------------------%
c
            call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
            call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_SBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c 
         else if ( type .eq. 5 ) then
c
c           %-------------------------------% 
c           | Perform y <-- OP*x            |
c           |    = inv[A-SIGMA*M]*A*x       |
c           | B*x = A*x has been computed   |
c           | and saved in workd(ipntr(3)). |
c           %-------------------------------%
c
            call dcopy (n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _SBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 6) then
c
c           %---------------------------------%
c           | Perform y <-- OP*x              |
c           | = inv[A-SIGMA*M]*(A+SIGMA*M)*x. | 
c           | (M*x) has been saved in         |
c           | workd(ipntr(3)).                |
c           %---------------------------------%
c
            if ( bmat .eq. 'G' ) then
               call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                    ab(itop,1), lda, workd(ipntr(1)), 1, 
     &                    zero, workd(ipntr(2)), 1)
               call daxpy( n, sigma, workd(ipntr(3)), 1, 
     &                    workd(ipntr(2)), 1 )
            else 
               call dcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call dgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                    lda, workd(ipntr(1)), 1, sigma, 
     &                    workd(ipntr(2)), 1)
            end if
            call dgbtrs('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
c
         end if
c
      else if (ido .eq. 2) then
c
c        %----------------------------------%
c        |        Perform y <-- B*x         | 
c        | Note when Buckling mode is used, |
c        | B = A, otherwise B=M.            | 
c        %----------------------------------%
c
         if (type .eq. 5) then
c
c           %---------------------%
c           | Buckling Mode, B=A. |
c           %---------------------%
c
            call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                ab(itop,1), lda, workd(ipntr(1)), 1, 
     &                zero, workd(ipntr(2)), 1)
         else

            call dgbmv('Notranspose', n, n, kl, ku, one, 
     &                mb(itop,1), lda, workd(ipntr(1)), 1, 
     &                zero, workd(ipntr(2)), 1)
         end if
c
      else 
c
c        %-----------------------------------------%
c        | Either we have convergence, or there is | 
c        | error.                                  |
c        %-----------------------------------------%
c
         if ( info .lt. 0) then
c
c           %--------------------------%
c           | Error message, check the |
c           | documentation in DSAUPD  |
c           %--------------------------%
c
            print *, ' '
            print *, ' Error with _saupd info = ',info
            print *, ' Check the documentation of _saupd '
            print *, ' '
            go to 9000
c
         else 
c
            if ( info .eq. 1) then
               print *, ' '
               print *, ' Maximum number of iterations reached.'
               print *, ' '
            else if ( info .eq. 3) then
               print *, ' '
               print *, ' No shifts could be applied during implicit',
     &                  ' Arnoldi update, try increasing NCV.'
               print *, ' '
            end if
c
            if (iparam(5) .gt. 0) then
c
               call dseupd ( rvec, 'A', select, d, z, ldz, sigma, 
     &                  bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     &                  iparam, ipntr, workd, workl, lworkl, info )            
c
               if ( info .ne. 0) then
c 
c                 %------------------------------------%
c                 | Check the documentation of dneupd. |
c                 %------------------------------------%
c
                  print *, ' ' 
                  print *, ' Error with _neupd = ', info
                  print *, ' Check the documentation of _neupd '
                  print *, ' ' 
                  go to 9000
c 
               end if
c
            end if
c
         end if
c
         go to 9000
c
      end if
c
c     %----------------------------------------%
c     | L O O P  B A C K to call DSAUPD again. |
c     %----------------------------------------%
c
      go to 90 
c
 9000 continue
c
      end
