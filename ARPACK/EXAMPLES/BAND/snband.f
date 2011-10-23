c \BeginDoc
c
c \Name: snband
c
c \Description:
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  Matrices A and B are stored in LAPACK-style banded form.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  snband can be called with one of the following modes:
c
c  Mode 1:  A*z = lambda*z.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*z = lambda*M*z, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c
c  Mode 3:  A*z = lambda*M*z, M symmetric semi-definite
c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*z = amu*z, then 
c           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
c           Note: If sigma is real, i.e. imaginary part of sigma is zero;
c                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M 
c                 amu == 1/(lambda-sigma). 
c  
c  Mode 4:  A*z = lambda*M*z, M symmetric semi-definite
c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*z = amu*z, then 
c           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
c
c
c  The choice of mode must be specified in IPARAM(7) defined below.
c
c \Usage
c   call snband
c      ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, 
c        WORKEV, V, N, AB, MB, LDA, RFAC, CFAC, KL, KU, WHICH, 
c        BMAT, NEV, TOL, RESID, NCV, V, LDV, IPARAM, WORKD, 
c        WORKL, LWORKL, WORKC, IWORK, INFO )
c
c \Arguments
c 
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether a basis for the invariant subspace corresponding 
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
c                                See Remarks below. 
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors; 
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c
c  DR      Real array of dimension NEV+1.  (OUTPUT)
c          On exit, DR contains the real part of the Ritz value approximations 
c          to the eigenvalues of A*z = lambda*B*z. 
c
c  DI      Real array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value 
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR. 
c
c          NOTE: When Ritz values are complex, they will come in complex 
c                conjugate pairs.  If eigenvectors are requested, the 
c                corresponding Ritz vectors will also come in conjugate 
c                pairs and the real and imaginary parts of these are 
c                represented in two consecutive columns of the array Z 
c                (see below).
c
c  Z       Real N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
c          On exit,
c          if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
c          Z represent approximate eigenvectors (Ritz vectors) corresponding
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z computed by SNAUPD.
c
c          The complex Ritz vector associated with the Ritz value
c          with positive imaginary part is stored in two consecutive
c          columns.  The first column holds the real part of the Ritz
c          vector and the second column holds the imaginary part.  The
c          Ritz vector associated with the Ritz value with negative
c          imaginary part is simply the complex conjugate of the Ritz vector
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by SNAUPD.  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT) 
c          The leading dimension of the array Z.  If Ritz vectors are 
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.  
c 
c  SIGMAR  Real  (INPUT) 
c          If IPARAM(7) = 3 or 4, represents the real part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.  
c 
c  SIGMAI  Real  (INPUT) 
c          If IPARAM(7) = 3 or 4, represents the imaginary part of the 
c          shift. 
c          Not referenced if IPARAM(7) = 1 or 2.  
c 
c  WORKEV  Real work array of dimension 3*NCV.  (WORKSPACE) 
c 
c  N       Integer.  (INPUT) 
c          Dimension of the eigenproblem.  
c 
c  AB      Real array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      Real array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set. 
c          The j-th column of M is stored in the j-th column of the
c          array AB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7) = 1
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, RFAC and CFAC. 
c
c  RFAC    Real array of LDA by N. (WORKSPACE/OUTPUT)
c          RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when IPARAM(7) = 3 is invoked with a real shift.
c          It is not referenced when IPARAM(7) = 1 or 4.
c
c  CFAC    Complex array of LDA by N. (WORKSPACE/OUTPUT)
c          CFAC is used to store (A-SIGMA*M) and its LU factors
c          when IPARAM(7) = 3 or 4 are used with a complex shift SIGMA.  
c          On exit, it contains the LU factors of (A-SIGMA*M).  
c          It is not referenced when IPARAM(7) = 1 or 2.
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
c            'LR' -> want the NEV eigenvalues of largest real part.
c            'SR' -> want the NEV eigenvalues of smallest real part.
c            'LI' -> want the NEV eigenvalues of largest imaginary part.
c            'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c          When IPARAM(7) = 3 or 4, WHICH should be set to 'LM' only. 
c          
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*z = lambda*z
c          BMAT = 'G' -> generalized eigenvalue problem A*z = lambda*M*z

c  NEV     Integer. (INPUT)
c          Number of eigenvalues to be computed.
c   
c  TOL     Real scalar.  (INPUT)
c          Stopping criteria: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = SLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine SLAMCH).
c
c  RESID   Real array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          Represents the dimension of the Arnoldi basis constructed
c          by snaupd for OP.
c
c  V       Real array N by NCV+1.  (OUTPUT)
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
c                       represent approximate Schur vectors that span the 
c                       desired invariant subspace.
c          NOTE: The array Z may be set equal to first NEV+1 columns of the 
c          Arnoldi basis vector array V computed by SNAUPD. In this case
c          if RVEC = .TRUE. and HOWMNY='A', then the first NCONV=IPARAM(5) 
c          are the desired Ritz vectors.
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
c           ----------------------------------------------------------
c          ISHIFT = 1: exact shift with respect to the current
c                      Hessenberg matrix H.  This is equivalent to
c                      restarting the iteration from the beginning
c                      after updating the starting vector with a linear
c                      combination of Ritz vectors associated with the
c                      "wanted" eigenvalues.
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
c
c          IPARAM(6) = IUPD
c          Not referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = IPARAM(7): 
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4; See under \Description of snband for the 
c          four modes available.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*z operations,
c                  NUMOPB = total number of B*z operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c WORKD    Real work array of length at least 3*n. (WORKSPACE)
c
c WORKL    Real work array of length LWORKL. (WORKSPACE)
c
c LWORKL   Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 6*NCV.
c
c WORKC    Complex array of length N. (WORKSPACE)
c          Workspace used when IPARAM(7) = 3 or 4 for storing a temporary 
c          complex vector.
c
c IWORK    Integer array of dimension at least N. (WORKSPACE)
c          Used when IPARAM(7)=2,3,4 to store the pivot information in the 
c          factorization of M or (A-SIGMA*M).
c            
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: The Schur form computed by LAPACK routine slahqr
c                could not be reordered by LAPACK routine strsen.
c                Re-enter subroutine SNEUPD with IPARAM(5)=NCV and 
c                increase the size of the arrays DR and DI to have 
c                dimension at least NCV and allocate at least NCV 
c                columns for Z. NOTE: Not necessary if Z and V share 
c                the same space. Please notify the authors.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine slahqr.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine strevc.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P'
c          = -14: SNAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: Overflow occurs when we try to transform the Ritz 
c                 values returned from SNAUPD to those of the original
c                 problem using Rayleigh Quotient. 
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
c     snaupd  ARPACK reverse communication interface routine.
c     sneupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     sgbtrf  LAPACK band matrix factorization routine.
c     sgbtrs  LAPACK band linear system solve routine.
c     cgbtrf  LAPACK complex band matrix factorization routine.
c     cgbtrs  LAPACK complex linear system solve routine.
c     slacpy  LAPACK matrix copy routine.
c     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     slamch  LAPACK routine to compute the underflow threshold.
c     scopy   Level 1 BLAS that copies one vector to another.
c     sdot    Level 1 BLAS that computes the dot product of two vectors.
c     snrm2   Level 1 BLAS that computes the norm of a vector.
c     sgbmv   Level 2 BLAS that computes the band matrix vector product.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let X' denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the real 
c     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; 
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
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
c FILE: nband.F   SID: 2.3   DATE OF SID: 10/17/00   RELEASE: 2
c
c\EndLib
c
c---------------------------------------------------------------------
c
      subroutine snband( rvec, howmny, select, dr, di, z, ldz,  sigmar, 
     &           sigmai, workev, n, ab, mb, lda, rfac,  cfac, kl, ku, 
     &           which, bmat, nev, tol, resid,  ncv, v, ldv, 
     &           iparam, workd, workl, lworkl, workc, iwork, info)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c 
      character        which*2, bmat, howmny
      integer          n, lda, kl, ku, nev, ncv, ldv,
     &                 ldz, lworkl, info  
      Real
     &                 tol, sigmar, sigmai 
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer          iparam(*), iwork(*)
      logical          select(*)
      Real
     &                 dr(*), di(*), resid(*), v(ldv,*), z(ldz,*),
     &                 ab(lda,*), mb(lda,*), rfac(lda,*), 
     &                 workd(*), workl(*), workev(*)
      Complex
     &                 cfac(lda,*), workc(*)
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
      Real        
     &                 numr, denr, deni, dmdul, safmin 
      logical          rvec, first 
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real
     &                  one, zero
      parameter        (one = 1.0E+0, zero = 0.0E+0)
c
c
c     %-----------------------------%
c     | LAPACK & BLAS routines used |
c     %-----------------------------%
c
      Real
     &                 sdot, snrm2, slapy2, slamch
      external         sdot, scopy, sgbmv, cgbtrf, cgbtrs, sgbtrf, 
     &                 sgbtrs, snrm2, slapy2, slacpy, slamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      Intrinsic        real, aimag, cmplx
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %--------------------------------%
c     | safmin = safe minimum is such  |
c     | that 1/sfmin does not overflow |
c     %--------------------------------%
c
      safmin = slamch('safmin')
c     
c     %----------------------------------------------------------------%
c     | Set type of the problem to be solved. Check consistency        |
c     | between BMAT and IPARAM(7).                                    |
c     | type = 1 --> Solving standard problem in regular mode.         |
c     | type = 2 --> Solving standard problem in shift-invert mode.    | 
c     | type = 3 --> Solving generalized problem in regular mode.      |
c     | type = 4 --> Solving generalized problem in shift-invert mode. |
c     | type = 5 --> Solving standard problem in shift-invert mode     |
c     |              using iparam(7) = 4 in SNAUPD.                    |
c     | type = 6 --> Solving generalized problem in shift-invert mode. | 
c     |              using iparam(7) = 4 in SNAUPD.                    |
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
      else if ( iparam(7) .eq. 4 .and. bmat .eq. 'I') then
         type = 5
      else if ( iparam(7) .eq. 4 .and. bmat .eq. 'G') then 
         type = 6
      else
         print*, ' '
         print*, 'BMAT is inconsistent with IPARAM(7).'
         print*, ' ' 
         go to 9000
      end if
c
c     %----------------------------------%       
c     | When type = 5,6 are used, sigmai |
c     | must be nonzero.                 |
c     %----------------------------------%
c
      if ( type .eq. 5 .or. type .eq. 6 ) then
          if ( sigmai .eq. zero ) then
             print*, ' '
             print*, '_NBAND: sigmai must be nonzero when type 5 or 6 
     &                is used. '
             print*, ' '
             go to 9000
          end if    
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
      if ( type .eq. 2 .or. type .eq. 5 ) then
c
c         %-------------------------------%
c         | Solving a standard eigenvalue |
c         | problem in shift-invert mode. |
c         | Factor (A-sigma*I).           |
c         %-------------------------------%
c
          if (sigmai .eq. zero) then
c            
c            %-----------------------------------%
c            | Construct (A-sigmar*I) and factor |
c            | in real arithmetic.               |
c            %-----------------------------------%
c
             call slacpy ('A', ibot, n, ab, lda, rfac, lda )
             do 10 j = 1, n
                rfac(imid,j) =  ab(imid,j) - sigmar
  10         continue
             call sgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr )
             if (ierr .ne. 0) then
                print*, ' ' 
                print*, ' _NBAND: Error with _gbtrf. '
                print*, ' '
                go to  9000
             end if
c
          else
c
c            %-----------------------------------%
c            | Construct (A-sigmar*I) and factor |
c            | in COMPLEX arithmetic.            |
c            %-----------------------------------%
c
             do 30 j = 1, n
                do 20 i = itop, ibot
                   cfac(i,j) = cmplx(ab(i,j))
  20            continue 
  30         continue
c
             do 40 j = 1, n
                cfac(imid,j) = cfac(imid,j) 
     $                         - cmplx(sigmar, sigmai)
  40         continue 
c 
             call cgbtrf(n, n, kl, ku, cfac, lda, iwork, ierr ) 
             if ( ierr .ne. 0) then
                print*, ' '
                print*, ' _NBAND: Error with _gbtrf. '
                print*, ' '
                go to  9000
             end if
c        
          end if
      
      else if ( type .eq. 3 ) then
c
c        %-----------------------------------------------%
c        | Solving generalized eigenvalue problem in     |
c        | regular mode. Copy M to rfac, and call LAPACK |
c        | routine sgbtrf to factor M.                   |
c        %-----------------------------------------------%
c
         call slacpy ('A', ibot, n, mb, lda, rfac, lda )
         call sgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
         if (ierr .ne. 0) then
             print*, ' ' 
             print*,'_NBAND:  Error with _gbtrf.'
             print*, ' '
             go to 9000
         end if
c
      else if ( type .eq. 4 .or. type .eq. 6 ) then
c
c        %-------------------------------------------%
c        | Solving generalized eigenvalue problem in |
c        | shift-invert mode.                        |
c        %-------------------------------------------%
c 
         if ( sigmai .eq. zero ) then
c
c            %--------------------------------------------%
c            | Construct (A - sigma*M) and factor in real |
c            | arithmetic.                                |
c            %--------------------------------------------%
c
             do 60 j = 1,n
                do 50 i = itop, ibot 
                   rfac(i,j) = ab(i,j) - sigmar*mb(i,j)
  50            continue
  60         continue
c
             call sgbtrf(n, n, kl, ku, rfac, lda, iwork, ierr)
             if ( ierr .ne. 0 )  then
                 print*, ' '
                 print*, '_NBAND: Error with _gbtrf.'
                 print*, ' '
                 go to 9000
             end if
c
         else
c
c            %-----------------------------------------------%
c            | Construct (A - sigma*M) and factor in complex |
c            | arithmetic.                                   |
c            %-----------------------------------------------% 
c
             do 80 j = 1,n
                do 70 i = itop, ibot 
                   cfac(i,j) = cmplx( ab(i,j)-sigmar*mb(i,j), 
     &                         -sigmai*mb(i,j) )
  70            continue 
  80         continue
c
             call cgbtrf(n, n, kl, ku, cfac, lda, iwork, ierr)
             if ( ierr .NE. 0 )  then
                print*, ' '
                print*, '_NBAND: Error with _gbtrf.'
                print*, ' '
                go to 9000
             end if
c 
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
      call snaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
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
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( type .eq. 2 ) then
c
            if (sigmai .eq. zero) then
c
c              %----------------------------------%
c              | Shift is real.  Perform          | 
c              | y <--- OP*x = inv[A-sigmar*I]*x  |
c              | to force the starting vector     |
c              | into the range of OP.            |
c              %----------------------------------%
c
               call scopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                       iwork, workd(ipntr(2)), n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' ' 
                  print*, ' _NBAND: Error with _bgtrs. '
                  print*, ' '
                  go to 9000
               end if
c
            else
c
c              %--------------------------------------------%
c              | Shift is COMPLEX. Perform                  |
c              | y <--- OP*x = Real_Part{inv[A-sigma*I]*x}  |
c              | to force the starting vector into the      | 
c              | range of OP.                               | 
c              %--------------------------------------------%
c
               do 100 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  100          continue
c
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                       iwork, workc, n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' '
                  print*, ' _NBAND: Error with _gbtrs. '
                  print*, ' '
                  go to 9000
               end if
c
               do 110 j = 1, n
                  workd(ipntr(2)+j-1) = real(workc(j))
  110          continue
c
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
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
c
            call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_NBAND: Error with _bgtrs.'
               print*, ' '
               go to 9000
            end if
c
         else if ( type .eq. 4 ) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x                      |
c           |         = Real_part{inv[A-SIGMA*M]*M}*x | 
c           | to force the starting vector into the   |
c           | range of OP.                            |
c           %-----------------------------------------%
c
            call sgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
            if ( sigmai .eq. zero ) then
c
c              %---------------------%
c              | Shift is real, stay |
c              | in real arithmetic. |
c              %---------------------%            
c
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                      iwork, workd(ipntr(2)), n, ierr)
               if (ierr .ne. 0) then
                  print*, ' ' 
                  print*, '_NBAND: Error with _gbtrs.'
                  print*, ' ' 
                  go to 9000
               end if
c
            else
c
c              %--------------------------%
c              | Goto complex arithmetic. |
c              %--------------------------%
c
               do 120 i = 1,n
                  workc(i) = cmplx(workd(ipntr(2)+i-1))
  120           continue 
c
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda, 
     &                      iwork, workc, n, ierr)
               if (ierr .ne. 0) then 
                  print*, ' '
                  print*, '_NBAND: Error with _gbtrs.' 
                  print*, ' '
                  go to 9000
               end if
c
               do  130 i = 1, n
                  workd(ipntr(2)+i-1) = real(workc(i))
  130          continue 
c
            end if
c
         else if ( type .eq. 5) then
c
c           %---------------------------------------% 
c           | Perform y <-- OP*x                    |
c           |    = Imaginary_part{inv[A-SIGMA*I]}*x |
c           | to force the starting vector into the |
c           | range of OP.                          |
c           %---------------------------------------%
c
            do 140 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  140       continue
c
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                    iwork, workc, n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _NBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
c
            do 150 j = 1, n
               workd(ipntr(2)+j-1) = aimag(workc(j))
  150       continue
c
         else if ( type .eq. 6 ) then
c
c           %----------------------------------------%
c           | Perform y <-- OP*x                     |
c           |       Imaginary_part{inv[A-SIGMA*M]*M} | 
c           | to force the starting vector into the  |
c           | range of OP.                           |
c           %----------------------------------------%
c
            call sgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
            do 160 i = 1,n
               workc(i) = cmplx(workd(ipntr(2)+i-1))
  160       continue 
c
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda, 
     &                   iwork, workc, n, ierr)
            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_NBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c
            do  170 i = 1, n
               workd(ipntr(2)+i-1) = aimag(workc(i))
  170       continue 
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
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( type .eq. 2) then
c
            if ( sigmai .eq. zero) then
c
c              %----------------------------------%
c              | Shift is real.  Perform          |
c              | y <--- OP*x = inv[A-sigmar*I]*x. |
c              %----------------------------------%
c
               call scopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                       iwork, workd(ipntr(2)), n, ierr)
            else
c
c              %------------------------------------------%
c              | Shift is COMPLEX. Perform                |
c              | y <-- OP*x = Real_Part{inv[A-sigma*I]*x} |
c              | in COMPLEX arithmetic.                   |
c              %------------------------------------------%
c
               do 180 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  180          continue
c
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                       iwork, workc, n, ierr)
               if ( ierr .ne. 0 ) then
                  print*, ' '
                  print*, '_NBAND: Error with _gbtrs.' 
                  print*, ' '
                  go to 9000
               end if
c
               do 190 j = 1, n
                  workd(ipntr(2)+j-1) = real(workc(j))
  190          continue
c
            end if
c
         else if ( type .eq. 3 ) then
c
c           %-----------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x |
c           %-----------------------------------%
c
            call sgbmv('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
c
            call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_NBAND: Error with _bgtrs.'
               print*, ' ' 
               go to 9000
            end if
c
         else if ( type .eq. 4 ) then
c
c           %--------------------------------------%
c           | Perform  y <-- inv(A-sigma*M)*(M*x). |
c           | (M*x) has been computed and stored   |
c           | in workd(ipntr(3)).                  |           
c           %--------------------------------------%
c
            if ( sigmai .eq. zero ) then
c
c              %------------------------%
c              | Shift is real, stay in |
c              | real arithmetic.       |
c              %------------------------%
c
               call scopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
               call sgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda, 
     &                       iwork, workd(ipntr(2)), n, ierr)
               if (ierr .ne. 0) then 
                  print*, ' '
                  print*, '_NBAND: Error with _gbtrs.' 
                  print*, ' '
                  go to 9000
               end if
c 
            else 
c
c              %---------------------------%
c              | Go to COMPLEX arithmetic. |
c              %---------------------------%
c
               do 200 i = 1,n
                  workc(i) = cmplx(workd(ipntr(3)+i-1))
  200          continue 
c
               call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda, 
     &                       iwork, workc, n, ierr)
               if (ierr .ne. 0) then 
                  print*, ' '
                  print*, '_NBAND: Error in _gbtrs.' 
                  print*, ' ' 
                  go to 9000
               end if
c
               do 210 i = 1,n
                  workd(ipntr(2)+i-1) = real(workc(i))
  210          continue 
c
            end if
c
         else if ( type .eq. 5 ) then
c
c           %---------------------------------------% 
c           | Perform y <-- OP*x                    |
c           |    = Imaginary_part{inv[A-SIGMA*I]*x} |
c           %---------------------------------------%
c
            do 220 j = 1, n
                  workc(j) = cmplx(workd(ipntr(1)+j-1))
  220       continue
c
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda,
     &                    iwork, workc, n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' _NBAND: Error with _gbtrs. '
               print*, ' '
               go to 9000
            end if
c
            do 230 j = 1, n
               workd(ipntr(2)+j-1) = aimag(workc(j))
  230       continue
c
         else if ( type .eq. 6) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x                      |
c           |   = Imaginary_part{inv[A-SIGMA*M]*M}*x. | 
c           %-----------------------------------------%
c
            do 240 i = 1,n
               workc(i) = cmplx(workd(ipntr(3)+i-1))
  240       continue 
c
            call cgbtrs ('Notranspose', n, kl, ku, 1, cfac, lda, 
     &                   iwork, workc, n, ierr)
            if (ierr .ne. 0) then 
               print*, ' '
               print*, '_NBAND: Error with _gbtrs.' 
               print*, ' '
               go to 9000
            end if
c
            do  250 i = 1, n
               workd(ipntr(2)+i-1) = aimag(workc(i))
  250       continue 
c
         end if
c
      else if (ido .eq. 2) then
c
c        %--------------------%
c        | Perform y <-- M*x  |
c        | Not used when      |
c        | type = 1,2.        |
c        %--------------------%
c
          call sgbmv('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &                lda, workd(ipntr(1)), 1, zero, 
     &                workd(ipntr(2)), 1)
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
c           | documentation in SNAUPD  |
c           %--------------------------%
c
            print *, ' '
            print *, ' Error with _naupd info = ',info
            print *, ' Check the documentation of _naupd '
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
               call sneupd ( rvec, 'A', select, dr, di, z, ldz, 
     &                 sigmar, sigmai, workev, bmat, n, which, 
     &                 nev, tol, resid, ncv, v, ldv, iparam,
     &                 ipntr, workd, workl, lworkl, info )            
c
               if ( info .ne. 0) then
c 
c                 %------------------------------------%
c                 | Check the documentation of SNEUPD. |
c                 %------------------------------------%
c
                  print *, ' ' 
                  print *, ' Error with _neupd = ', info
                  print *, ' Check the documentation of _neupd '
                  print *, ' ' 
                  go to 9000
c
               else if ( sigmai .ne. zero ) then 
c 
                  if ( type .eq. 4 .or. type .eq. 6 ) then
c 
                     first = .true.
                     do 270 j = 1, iparam(5) 
c
c                    %----------------------------------%
c                    | Use Rayleigh Quotient to recover |
c                    | eigenvalues of the original      |
c                    | generalized eigenvalue problem.  |
c                    %----------------------------------%
c
                     if ( di(j) .eq. zero ) then
c
c                       %--------------------------------------%
c                       | Eigenvalue is real. Compute          |
c                       | d = (x'*inv[A-sigma*M]*M*x) / (x'*x) | 
c                       %--------------------------------------%
c
                        call sgbmv('Nontranspose', n, n, kl, ku, one, 
     $                     mb(itop,1), lda, z(1,j), 1, zero,
     $                     workd, 1)
                        do i = 1, n
                           workc(i) = cmplx(workd(i))
                        end do
                        call cgbtrs ('Notranspose', n, kl, ku, 1, 
     $                        cfac, lda, iwork, workc, n, info)
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n, z(1,j), 1, workd, 1)
                        deni = sdot(n, z(1,j), 1, workd(n+1), 1)
                        numr  = snrm2(n, z(1,j), 1)**2
                        dmdul = slapy2(denr,deni)**2
                        if ( dmdul .ge. safmin ) then
                           dr(j) = sigmar + numr*denr / dmdul
                        else
c
c                          %---------------------%
c                          | dmdul is too small. |
c                          | Exit to avoid       |
c                          | overflow.           |
c                          %---------------------%
c
                           info = -15
                           go to 9000
                        end if
c
                     else if (first) then 
c
c                       %------------------------%
c                       | Eigenvalue is complex. |
c                       | Compute the first one  |
c                       | of the conjugate pair. |
c                       %------------------------%
c
c                       %-------------%
c                       | Compute M*x |
c                       %-------------%
c
                        call sgbmv('Nontranspose', n, n, kl, ku,
     $                      one, mb(itop,1), lda, z(1,j), 1, zero,
     $                      workd, 1)
                        call sgbmv('Nontranspose', n, n, kl, ku, 
     $                       one, mb(itop,1), lda, z(1,j+1), 1,
     $                       zero, workd(n+1), 1)
                        do i = 1, n
                           workc(i) = cmplx(workd(i),workd(i+n))
                        end do
c
c                       %----------------------------%
c                       | Compute inv(A-sigma*M)*M*x |
c                       %----------------------------%
c
                        call cgbtrs('Notranspose',n,kl,ku,1,cfac, 
     $                     lda, iwork, workc, n, info)
c
c                       %-------------------------------%
c                       | Compute x'*inv(A-sigma*M)*M*x |
c                       %-------------------------------%
c
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n,z(1,j),1,workd,1)
                        denr = denr+sdot(n,z(1,j+1),1,workd(n+1),1)
                        deni = sdot(n,z(1,j),1,workd(n+1),1)
                        deni = deni - sdot(n,z(1,j+1),1,workd,1)
c
c                       %----------------%
c                       | Compute (x'*x) |
c                       %----------------%
c
                        numr = slapy2( snrm2(n, z(1,j), 1), 
     &                         snrm2(n, z(1, j+1), 1) )**2
c
c                       %----------------------------------------%
c                       | Compute (x'x) / (x'*inv(A-sigma*M)*Mx) |
c                       %----------------------------------------%
c
                        dmdul = slapy2(denr,deni)**2
                        if ( dmdul .ge. safmin ) then
                           dr(j) = sigmar+numr*denr / dmdul
                           di(j) = sigmai-numr*deni / dmdul
                           first = .false.
                        else
c
c                          %---------------------%
c                          | dmdul is too small. |
c                          | Exit to avoid       |
c                          | overflow.           |
c                          %---------------------%
c
                           info = -15
                           go to 9000
c
                        end if
c
                     else
c
c                       %---------------------------%
c                       | Get the second eigenvalue |
c                       | of the conjugate pair by  |
c                       | taking the conjugate of   |
c                       | previous one.             |
c                       %---------------------------%
c
                        dr(j) = dr(j-1)
                        di(j) = -di(j-1)
                        first = .true.
c
                     end if
c
  270                continue 
c
                  else if ( type .eq. 2 .or. type .eq. 5) then
c
                     first = .true.
                     do 280 j = 1, iparam(5)
c
c                    %----------------------------------%
c                    | Use Rayleigh Quotient to recover |
c                    | eigenvalues of the original      |
c                    | standard eigenvalue problem.     |
c                    %----------------------------------%
c
                     if ( di(j) .eq. zero ) then
c
c                       %-------------------------------------%
c                       | Eigenvalue is real. Compute         |
c                       | d = (x'*inv[A-sigma*I]*x) / (x'*x). |
c                       %-------------------------------------%
c
                        do i = 1, n
                           workc(i) = cmplx(z(i,j))
                        end do
                        call cgbtrs ('Notranspose', n, kl, ku, 1, 
     $                        cfac, lda, iwork, workc, n, info)
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n,z(1,j),1,workd,1)
                        deni = sdot(n,z(1,j),1,workd(n+1),1)
                        numr  = snrm2(n, z(1,j), 1)**2
                        dmdul = slapy2(denr,deni)**2
                        if ( dmdul .ge. safmin ) then
                           dr(j) = sigmar + numr*denr / dmdul
                        else
c
c                          %---------------------%
c                          | dmdul is too small. |
c                          | Exit to avoid       |
c                          | overflow.           |
c                          %---------------------%
c
                           info = -15
                           go to 9000
c
                        end if
c
                     else if (first) then
c
c                       %------------------------%
c                       | Eigenvalue is complex. |
c                       | Compute the first one  |
c                       | of the conjugate pair. |
c                       %------------------------%
c
                        do i = 1, n
                           workc(i) = cmplx( z(i,j), z(i,j+1) )
                        end do
c
c                       %---------------------------%
c                       | Compute inv[A-sigma*I]*x. |
c                       %---------------------------%
c
                        call cgbtrs('Notranspose',n,kl,ku,1,cfac,
     $                       lda, iwork, workc, n, info)
c
c                       %-----------------------------%
c                       | Compute x'*inv(A-sigma*I)*x |
c                       %-----------------------------%
c
                        do i = 1, n
                           workd(i) = real(workc(i))
                           workd(i+n) = aimag(workc(i))
                        end do
                        denr = sdot(n,z(1,j),1,workd,1)
                        denr = denr+sdot(n,z(1,j+1),1,workd(n+1),1)
                        deni = sdot(n,z(1,j),1,workd(n+1),1)
                        deni = deni - sdot(n,z(1,j+1),1,workd,1)
c
c                       %----------------%
c                       | Compute (x'*x) |
c                       %----------------%
c
                        numr = slapy2( snrm2(n, z(1,j), 1),
     &                         snrm2(n, z(1,j+1), 1))**2
c
c                       %----------------------------------------%
c                       | Compute (x'x) / (x'*inv(A-sigma*I)*x). |
c                       %----------------------------------------%
c
                        dmdul = slapy2(denr,deni)**2
                        if (dmdul .ge. safmin) then  
                           dr(j) = sigmar+numr*denr / dmdul
                           di(j) = sigmai-numr*deni / dmdul
                           first = .false.
                        else
c
c                          %---------------------%
c                          | dmdul is too small. |
c                          | Exit to avoid       |
c                          | overflow.           |
c                          %---------------------%
c
                           info = -15
                           go to 9000
                        end if
c
                     else
c
c                       %---------------------------%
c                       | Get the second eigenvalue |
c                       | of the conjugate pair by  |
c                       | taking the conjugate of   |
c                       | previous one.             |
c                       %---------------------------%
c
                        dr(j) = dr(j-1)
                        di(j) = -di(j-1)
                        first = .true.
c
                     end if
c
  280                continue
c
                  end if
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
c     | L O O P  B A C K to call SNAUPD again. |
c     %----------------------------------------%
c
      go to 90 
c
 9000 continue
c
      end
