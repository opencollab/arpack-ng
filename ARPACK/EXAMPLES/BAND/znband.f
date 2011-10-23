c \BeginDoc
c
c \Name: znband 
c
c \Description:
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
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are commonly called Ritz values and Ritz vectors respectively.  They are 
c  referred to as such in the comments that follow.  The computed orthonormal 
c  basis for the invariant subspace corresponding to these Ritz values is 
c  referred to as a Schur basis. 
c
c  znband  can be called with one of the following modes:
c
c  Mode 1:  A*z = lambda*z.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*z = lambda*M*z, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c
c  Mode 3:  A*z = lambda*M*z, M symmetric semi-definite
c           ===> OP = inv[A - sigma*M]*M   and  B = M.
c           ===> shift-and-invert mode.
c
c  Choice of different modes can be specified in IPARAM(7) defined below.
c
c \Usage
c   call znband 
c      ( RVEC, HOWMNY, SELECT, D , Z, LDZ, SIGMA, WORKEV, N, AB, 
c        MB, LDA, FAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, 
c        V, LDV, IPARAM, WORKD, WORKL, LWORKL, RWORK, IWORK, INFO )
c
c \Arguments
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether a basis for the invariant subspace corresponding
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
c                                See Remarks below.
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies the form of the invariant subspace to be computed 
c          corresponding to the converged Ritz values.
c          = 'A': Compute NEV Ritz vectors;
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the real Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT need not be initialized
c          but it is used as internal workspace.
c
c  D       Complex*16  array of dimension NEV+1.  (OUTPUT)
c          On exit, D contains the  Ritz  approximations
c          to the eigenvalues lambda for A*z = lambda*B*z.
c
c  Z       Complex*16  N by NEV array     (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
c          Z represents approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z.
c
c          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required, 
c          the array Z may be set equal to first NEV columns of the 
c          array V.  
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ) is required.
c          In any case,  LDZ .ge. 1 is required.
c
c  SIGMA   Complex*16   (INPUT)
c          If IPARAM(7) = 3 then SIGMA represents the shift.
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  WORKEV  Complex*16  work array of dimension NCV.  (WORKSPACE)
c 
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  AB      Complex*16  array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      Complex*16  array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set. 
c          The j-th column of M is stored in the j-th column of the
c          array MB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7)=1.
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, FAC.
c
c  FAC     Complex*16  array of LDA by N. (WORKSPACE/OUTPUT)
c          FAC is used to store the LU factors of MB when mode 2
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when mode 3 is invoked.
c          It is not referenced when IPARAM(7)=1.
c
c  KL      Integer. (INPUT)
c          Max(number of subdiagonals of A, number of subdiagonals of M)
c
c  KU      Integer. (OUTPUT)
c          Max(number of superdiagonals of A, number of superdiagonals of M)
c
c  WHICH   Character*2.  (INPUT)
c          When mode 1,2 are used, WHICH can be set to any one of
c          the following.
c  
c            'LM' -> want the NEV eigenvalues of largest magnitude.
c            'SM' -> want the NEV eigenvalues of smallest magnitude.
c            'LR' -> want the NEV eigenvalues of largest real part.
c            'SR' -> want the NEV eigenvalues of smallest real part.
c            'LI' -> want the NEV eigenvalues of largest imaginary part.
c            'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c          When mode 3 is used, WHICH should be set to 'LM' only. 
c          
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

c  NEV     Integer. (INPUT)
c          Number of eigenvalues of to be computed.
c   
c  TOL     Double precision  scalar.  (INPUT)
c          Stopping criteria: the relative accuracy of the Ritz value
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = dlamch ('EPS')  (machine precision as computed
c                    by the LAPACK auxilliary subroutine dlamch ).
c
c  RESID   Complex*16  array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 2 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. 
c
c  V       Complex*16  array N by NCV.  (OUTPUT)
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then 
c          the first NCONV=IPARAM(5) columns of V will contain Ritz vectors 
c          of the eigensystem A*z = lambda*B*z.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.  LDV must be great than or equal to N.
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
c          IPARAM(2) = Not referenced.
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
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2 or 3; See under \Description of znband  for the 
c          three modes available.
c
c WORKD    Complex*16  work array of length at least 3*n. (WORKSPACE)
c
c WORKL    Complex*16  work array of length LWORKL. (WORKSPACE) 
c
c LWORKL   Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 5*NCV.
c
c RWORK    Double precision  array of length N (WORKSPACE)
c          Workspace used in znaupd .
c
c IWORK    Integer array of dimension at least N. (WORKSPACE)
c          Used to mode 2,3. Store the pivot information in the 
c          factorization of M or (A-SIGMA*M).
c            
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation.
c                This should never happened.
c          = -10: IPARAM(7) must be 1,2,3.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: ZNAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c
c \EndDoc
c
c------------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called
c     znaupd   ARPACK reverse communication interface routine.
c     zneupd   ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     zgbtrf   LAPACK band matrix factorization routine.
c     zgbtrs   LAPACK band linear system solve routine.
c     zlacpy   LAPACK matrix copy routine.
c     zcopy    Level 1 BLAS that copies one vector to another.
c     dznrm2   Level 1 BLAS that computes the norm of a vector.
c     zgbmv    Level 2 BLAS that computes the band matrix vector product.
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Ph.D thesis, TR95-13, Rice Univ,
c     May 1995.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
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
c-----------------------------------------------------------------------
c
      subroutine znband (rvec, howmny, select, d , z, ldz, sigma,
     &                 workev, n, ab, mb, lda, fac, kl, ku, which, 
     &                 bmat, nev, tol, resid, ncv, v, ldv, iparam, 
     &                 workd, workl, lworkl, rwork, iwork, info )
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c 
      Character        which*2, bmat, howmny
      Logical          rvec
      Integer          n, lda, kl, ku, nev, ncv, ldv,
     &                 ldz, lworkl, info  
      Complex*16          
     &                 sigma 
      Double precision 
     &                 tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Integer          iparam(*), iwork(*)
      Logical          select(*)
      Complex*16          
     &                 d(*), resid(*), v(ldv,*), z(ldz,*),
     &                 ab(lda,*), mb(lda,*), fac(lda,*), 
     &                 workd(*), workl(*), workev(*)
      Double precision 
     &                 rwork(*)
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
      integer          ido, i, j, mode, ierr, itop, imid, ibot
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16          
     &                  one, zero
      parameter        (one  = (1.0D+0, 0.0D+0) ,
     &                  zero = (0.0D+0, 0.0D+0) )
c
c     %-----------------------------%
c     | LAPACK & BLAS routines used |
c     %-----------------------------%
c
      Double precision 
     &                 dznrm2 
      external         zcopy , zgbmv , zgbtrf , zgbtrs , dznrm2 , zlacpy 
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c     
      mode = iparam(7)
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
      if ( mode .eq. 2 ) then
c
c         %-----------------------------------------------%
c         | Copy M to fac and Call LAPACK routine zgbtrf   |
c         | to factor M.                                  |
c         %-----------------------------------------------%
c
          call zlacpy  ('A', ibot, n, mb, lda, fac, lda )
          call zgbtrf (n, n, kl, ku, fac, lda, iwork, ierr) 
          if (ierr .ne. 0) then
              print*, ' ' 
              print*,'_band:  error in _gbtrf'
              print*, ' '
              go to 9000
          end if
c
      else if ( mode .eq. 3 ) then
c
          if (bmat .eq. 'I') then
c
c            %-------------------------%
c            | Construct (A - sigma*I) |
c            %-------------------------%
c
             call zlacpy  ('A', ibot, n, ab, lda, fac, lda )
             do 10 j = 1,n
                fac(imid,j) = ab(imid,j) - sigma
  10         continue
c
          else
c
c            %---------------------------%
c            | Construct (A - sigma*M)   |
c            %---------------------------%
c
             do 30 j = 1,n
                do 20 i = itop, ibot 
                   fac(i,j) = ab(i,j) - sigma*mb(i,j)
  20            continue
  30         continue
c
          end if
c
c         %------------------------%
c         | Factor (A - sigma*M)   |
c         %------------------------%
c
          call zgbtrf (n, n, kl, ku, fac, lda, iwork, ierr)
          if ( ierr .ne. 0 )  then
              print*, ' '
              print*, '_band: error in _gbtrf.'
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
  40  continue 
c
      call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv,
     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
     &              rwork,info )

c
      if (ido .eq. -1) then
c
         if ( mode .eq. 1) then
c
c           %----------------------------%
c           | Perform  y <--- OP*x = A*x |
c           %----------------------------%
c
            call zgbmv ('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( mode .eq. 2 ) then
c
c           %-----------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x |
c           %-----------------------------------%
c
            call zgbmv ('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
c
            call zgbtrs  ('Notranspose', n, kl, ku, 1, fac, lda, 
     &                    iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_band: error in sbgtrs.'
               print*, ' '
               go to 9000
            end if
c
         else if ( mode .eq. 3 ) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x                      |
c           |           = inv[A-SIGMA*M]*M* x 
c           | to force the starting vector into the   |
c           | range of OP.                            |
c           %-----------------------------------------%
c
            call zgbmv ('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
            call zgbtrs  ('Notranspose', n, kl, ku, 1, fac, lda, 
     &                   iwork, workd(ipntr(2)), n, ierr)
            if (ierr .ne. 0) then
               print*, ' ' 
               print*, '_band: error in _gbtrs.'
               print*, ' ' 
               go to 9000
            end if
c
         end if
c
      else if (ido .eq. 1) then
c
         if ( mode .eq. 1) then
c
c           %----------------------------%
c           | Perform  y <--- OP*x = A*x |
c           %----------------------------%
c
            call zgbmv ('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                 lda, workd(ipntr(1)), 1, zero, 
     &                 workd(ipntr(2)), 1)
c
         else if ( mode .eq. 2 ) then
c
c           %-----------------------------------%
c           | Perform  y <--- OP*x = inv[M]*A*x |
c           %-----------------------------------%
c
            call zgbmv ('Notranspose', n, n, kl, ku, one, ab(itop,1), 
     &                  lda, workd(ipntr(1)), 1, zero, 
     &                  workd(ipntr(2)), 1)
c
            call zgbtrs  ('Notranspose', n, kl, ku, 1, fac, lda, 
     &                    iwork, workd(ipntr(2)), ldv, ierr)
            if (ierr .ne. 0) then
               print*, ' '
               print*, '_band: error in sbgtrs.'
               print*, ' ' 
               go to 9000
            end if
c
         else if ( mode .eq. 3 ) then
c
            if ( bmat .eq. 'I' ) then
c
c              %----------------------------------%
c              | Perform  y <-- inv(A-sigma*I)*x. |
c              %----------------------------------%
c
               call zcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)
               call zgbtrs  ('Notranspose', n, kl, ku, 1, fac, lda,
     &                    iwork, workd(ipntr(2)), n, ierr)
               if (ierr .ne. 0) then
                  print*, ' '
                  print*, '_band: error in _gbtrs.'
                  print*, ' '
                  go to 9000
               end if
c
            else
c  
c              %--------------------------------------%
c              | Perform  y <-- inv(A-sigma*M)*(M*x). |
c              | (M*x) has been computed and stored   |
c              | in workd(ipntr(3)).                  |           
c              %--------------------------------------%
c
               call zcopy (n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
               call zgbtrs  ('Notranspose', n, kl, ku, 1, fac, lda, 
     &                      iwork, workd(ipntr(2)), n, ierr)
               if (ierr .ne. 0) then 
                  print*, ' '
                  print*, '_band: error in _gbtrs.' 
                  print*, ' '
                  go to 9000
               end if
c
            end if
c
         endif
c
      else if (ido .eq. 2) then
c
c        %--------------------%
c        | Perform y <-- M*x  |
c        %--------------------%
c
          call zgbmv ('Notranspose', n, n, kl, ku, one, mb(itop,1), 
     &                lda, workd(ipntr(1)), 1, zero, 
     &                workd(ipntr(2)), 1)
c
      else 
c
c        %-------------------------------------------%
c        |   Either we have convergence, or there is | 
c        |   error.                                  |
c        %-------------------------------------------%
c
         if ( info .ne. 0) then
c
c           %--------------------------%
c           | Error message, check the |
c           | documentation in dnaupd  |
c           %--------------------------%
c
            print *, ' '
            print *, ' Error with _naupd info = ',info
            print *, ' Check the documentation of _naupd '
            print *, ' '
c
         else 
c
            call zneupd  (rvec, howmny , select, d, z, ldz, sigma,
     &                   workev, bmat, n, which, nev, tol,
     &                   resid, ncv, v, ldv, iparam, ipntr, workd,
     &                   workl, lworkl, rwork, info)
c
            if ( info .ne. 0) then
c 
c              %------------------------------------%
c              | Check the documentation of zneupd . |
c              %------------------------------------%
c
               print *, ' ' 
               print *, ' Error with _neupd = ', info
               print *, ' Check the documentation of _neupd '
               print *, ' ' 
c 
            endif 
c
         end if
c
         go to 9000
c
      end if
c
c     %----------------------------------------%
c     | L O O P  B A C K to call znaupd  again. |
c     %----------------------------------------%
c
      go to 40 
c
 9000 continue
c
      end
