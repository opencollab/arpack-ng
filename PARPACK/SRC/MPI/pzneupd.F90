!\BeginDoc
!
!\Name: pzneupd
!
! Message Passing Layer: MPI
!
!\Description:
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  basis is always computed.  There is an additional storage cost of n*nev
!  if both are requested (in this case a separate array Z must be supplied).
!
!  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
!  are derived from approximate eigenvalues and eigenvectors of
!  of the linear operator OP prescribed by the MODE selection in the
!  call to PZNAUPD.  PZNAUPD must be called before this routine is called.
!  These approximate eigenvalues and vectors are commonly called Ritz
!  values and Ritz vectors respectively.  They are referred to as such
!  in the comments that follow.   The computed orthonormal basis for the
!  invariant subspace corresponding to these Ritz values is referred to as a
!  Schur basis.
!
!  The definition of OP as well as other terms and the relation of computed
!  Ritz values and vectors of OP with respect to the given problem
!  A*z = lambda*B*z may be found in the header of PZNAUPD.  For a brief
!  description, see definitions of IPARAM(7), MODE and WHICH in the
!  documentation of PZNAUPD.
!
!\Usage:
!  call pzneupd
!     ( COMM, RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT,
!       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD,
!       WORKL, LWORKL, RWORK, INFO )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether a basis for the invariant subspace corresponding
!          to the converged Ritz value approximations for the eigenproblem
!          A*z = lambda*B*z is computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
!                                See Remarks below.
!
!  HOWMNY  Character*1  (INPUT)
!          Specifies the form of the basis for the invariant subspace
!          corresponding to the converged Ritz values that is to be computed.
!
!          = 'A': Compute NEV Ritz vectors;
!          = 'P': Compute NEV Schur vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the  Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A' or 'P', SELECT need not be initialized
!          but it is used as internal workspace.
!
!  D       Complex*16 array of dimension NEV+1.  (OUTPUT)
!          On exit, D contains the  Ritz  approximations
!          to the eigenvalues lambda for A*z = lambda*B*z.
!
!  Z       Complex*16 N by NEV array    (OUTPUT)
!          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
!          Z represents approximate eigenvectors (Ritz vectors) corresponding
!          to the NCONV=IPARAM(5) Ritz values for eigensystem
!          A*z = lambda*B*z.
!
!          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
!
!          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
!          the array Z may be set equal to first NEV+1 columns of the Arnoldi
!          basis array V computed by PZNAUPD.  In this case the Arnoldi basis
!          will be destroyed and overwritten with the eigenvector basis.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ .ge.  max( 1, N ) is required.
!          In any case,  LDZ .ge. 1 is required.
!
!  SIGMA   Complex*16  (INPUT)
!          If IPARAM(7) = 3 then SIGMA represents the shift.
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  WORKEV  Complex*16 work array of dimension 2*NCV.  (WORKSPACE)
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to PZNAUPD that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, RWORK, INFO
!
!         must be passed directly to ZNEUPD following the last call
!         to PZNAUPD.  These arguments MUST NOT BE MODIFIED between
!         the the last call to PZNAUPD and the call to ZNEUPD.
!
!  Three of these parameters (V, WORKL and INFO) are also output parameters:
!
!  V       Complex*16 N by NCV array.  (INPUT/OUTPUT)
!
!          Upon INPUT: the NCV columns of V contain the Arnoldi basis
!                      vectors for OP as constructed by PZNAUPD .
!
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
!                       contain approximate Schur vectors that span the
!                       desired invariant subspace.
!
!          NOTE: If the array Z has been set equal to first NEV+1 columns
!          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
!          Arnoldi basis held by V has been overwritten by the desired
!          Ritz vectors.  If a separate array Z has been passed then
!          the first NCONV=IPARAM(5) columns of V will contain approximate
!          Schur vectors that span the desired invariant subspace.
!
!  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
!          PZNAUPD.  They are not changed by PZNEUPD.
!          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
!          untransformed Ritz values, the untransformed error estimates of
!          the Ritz values, the upper triangular matrix for H, and the
!          associated matrix representation of the invariant subspace for H.
!
!          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
!          of the above information computed by PZNEUPD.
!          -------------------------------------------------------------
!          IPNTR(9):  pointer to the NCV RITZ values of the
!                     original system.
!          IPNTR(10): Not used
!          IPNTR(11): pointer to the NCV corresponding error estimates.
!          IPNTR(12): pointer to the NCV by NCV upper triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     PZNEUPD if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!
!          =  1: The Schur form computed by LAPACK routine csheqr
!                could not be reordered by LAPACK routine ztrsen.
!                Re-enter subroutine pzneupd with IPARAM(5)=NCV and
!                increase the size of the array D to have
!                dimension at least dimension NCV and allocate at least NCV
!                columns for Z. NOTE: Not necessary if Z and V share
!                the same space. Please notify the authors if this error
!                occurs.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from LAPACK eigenvalue calculation.
!                This should never happened.
!          = -9: Error return from calculation of eigenvectors.
!                Informational error from LAPACK routine ztrevc.
!          = -10: IPARAM(7) must be 1,2,3
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: HOWMNY = 'S' not yet implemented
!          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
!          = -14: PZNAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: ZNEUPD got a different count of the number of converged
!                 Ritz values than ZNAUPD got.  This indicates the user
!                 probably made an error in passing data from ZNAUPD to
!                 ZNEUPD or that the data was modified before entering
!                 ZNEUPD.
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
!     "How to Implement the Spectral Transformation", Math Comp.,
!     Vol. 48, No. 178, April, 1987 pp. 664-673.
!
!\Routines called:
!     pivout  Parallel ARPACK utility routine that prints integers.
!     pzmout  Parallel ARPACK utility routine that prints matrices
!     pzvout  Parallel ARPACK utility routine that prints vectors.
!     zgeqr2  LAPACK routine that computes the QR factorization of
!             a matrix.
!     zlacpy  LAPACK matrix copy routine.
!     zlahqr  LAPACK routine that computes the Schur form of a
!             upper Hessenberg matrix.
!     zlaset  LAPACK matrix initialization routine.
!     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper triangular form.
!     ztrsen  LAPACK routine that re-orders the Schur form.
!     zunm2r  LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     pdlamch10 ScaLAPACK routine that determines machine constants.
!     ztrmm   Level 3 BLAS matrix times an upper triangular matrix.
!     zgeru   Level 2 BLAS rank one update to a matrix.
!     zcopy   Level 1 BLAS that copies one vector to another .
!     zscal   Level 1 BLAS that scales a vector.
!     zdscal  Level 1 BLAS that scales a complex vector by a real number.
!     dznrm2  Level 1 BLAS that computes the norm of a complex vector.
!
!\Remarks
!
!  1. Currently only HOWMNY = 'A' and 'P' are implemented.
!
!  2. Schur vectors are an orthogonal representation for the basis of
!     Ritz vectors. Thus, their numerical properties are often superior.
!     If RVEC = .true. then the relationship
!             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
!     V(:,1:IPARAM(5))` * V(:,1:IPARAM(5)) = I are approximately satisfied.
!     Here T is the leading submatrix of order IPARAM(5) of the
!     upper triangular matrix stored workl(ipntr(12)).
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Parallel Modifications
!     Kristi Maschhoff
!
!\Revision history:
!     Starting Point: Complex Serial Code FILE: neupd.F   SID: 2.2
!
!\SCCS Information:
! FILE: neupd.F   SID: 1.9   DATE OF SID: 10/25/03
!
!\EndLib
!
!-----------------------------------------------------------------------
      subroutine pzneupd
     &         ( comm , rvec  , howmny, select, d    ,
     &           z    , ldz   , sigma , workev, bmat ,
     &           n    , which , nev   , tol   , resid,
     &           ncv  , v     , ldv   , iparam, ipntr,
     &           workd, workl , lworkl, rwork , info )
!
!     %--------------------%
!     | MPI Communicator |
!     %--------------------%
!
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
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Complex*16
     &           sigma
      Double precision
     &           tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision
     &           rwork(ncv)
      Complex*16
     &           d(nev)     , resid(n)  , v(ldv,ncv)   ,
     &           z(ldz, nev), workd(3*n), workl(lworkl),
     &           workev(2*ncv)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex*16
     &           one, zero
      parameter  (one = (1.0, 0.0), zero = (0.0, 0.0))
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  type*6
      integer    bounds, ierr  , ih    , ihbds, iheig , nconv ,
     &           invsub, iuptri, iwev  , j    , ldh   , ldq   ,
     &           mode  , msglvl, ritz  , wr   , k     , irz   ,
     &           ibd   , outncv, iq    , np   , numcnv, jj    ,
     &           ishift
      Complex*16
     &           rnorm, temp, vl(1)
      Double precision
     &           conds, sep, rtemp, eps23
      logical    reord
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   zcopy ,zgeru,zgeqr2,zlacpy,pzmout,
     &           zunm2r,ztrmm,pzvout,pivout,
     &           zlahqr
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision
     &           dznrm2,pdlamch10,dlapy2
      external   dznrm2,pdlamch10,dlapy2
!
      Complex*16
     &           zzdotc
      external   zzdotc
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
      msglvl = mceupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
!
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
      eps23 = pdlamch10(comm, 'Epsilon-Machine')
      eps23 = eps23**(2.0 / 3.0)
!
!     %-------------------------------%
!     | Quick return                  |
!     | Check for incompatible input  |
!     %-------------------------------%
!
      ierr = 0
!
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1 .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 4*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
!
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 ) then
         type = 'SHIFTI'
      else
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
!
!     %--------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, WORKEV, Q   |
!     | etc... and the remaining workspace.                    |
!     | Also update pointer to be used on output.              |
!     | Memory is laid out as follows:                         |
!     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
!     | workl(ncv*ncv+1:ncv*ncv+ncv) := ritz values            |
!     | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := error bounds     |
!     %--------------------------------------------------------%
!
!     %-----------------------------------------------------------%
!     | The following is used and set by ZNEUPD.                  |
!     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := The untransformed |
!     |                                      Ritz values.         |
!     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
!     |                                      error bounds of      |
!     |                                      the Ritz values      |
!     | workl(ncv*ncv+4*ncv+1:2*ncv*ncv+4*ncv) := Holds the upper |
!     |                                      triangular matrix    |
!     |                                      for H.               |
!     | workl(2*ncv*ncv+4*ncv+1: 3*ncv*ncv+4*ncv) := Holds the    |
!     |                                      associated matrix    |
!     |                                      representation of    |
!     |                                      the invariant        |
!     |                                      subspace for H.      |
!     | GRAND total of NCV * ( 3 * NCV + 4 ) locations.           |
!     %-----------------------------------------------------------%
!
      ih     = ipntr(5)
      ritz   = ipntr(6)
      iq     = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheig  = bounds + ldh
      ihbds  = iheig  + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheig
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wr = 1
      iwev = wr + ncv
!
!     %-----------------------------------------%
!     | irz points to the Ritz values computed  |
!     |     by _neigh before exiting _naup2.    |
!     | ibd points to the Ritz estimates        |
!     |     computed by _neigh before exiting   |
!     |     _naup2.                             |
!     %-----------------------------------------%
!
      irz = ipntr(14) + ncv*ncv
      ibd = irz + ncv
!
!     %------------------------------------%
!     | RNORM is B-norm of the RESID(1:N). |
!     %------------------------------------%
!
      rnorm = workl(ih+2)
      workl(ih+2) = zero
!
      if (msglvl .gt. 2) then
         call pzvout(comm, logfil, ncv, workl(irz), ndigit,
     &   '_neupd: Ritz values passed in from _NAUPD.')
         call pzvout(comm, logfil, ncv, workl(ibd), ndigit,
     &   '_neupd: Ritz estimates passed in from _NAUPD.')
      end if
!
      if (rvec) then
!
         reord = .false.
!
!        %---------------------------------------------------%
!        | Use the temporary bounds array to store indices   |
!        | These will be used to mark the select array later |
!        %---------------------------------------------------%
!
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
!
!        %-------------------------------------%
!        | Select the wanted Ritz values.      |
!        | Sort the Ritz values so that the    |
!        | wanted ones appear at the tailing   |
!        | NEV positions of workl(irr) and     |
!        | workl(iri).  Move the corresponding |
!        | error estimates in workl(ibd)       |
!        | accordingly.                        |
!        %-------------------------------------%
!
         np     = ncv - nev
         ishift = 0
         call pzngets(comm, ishift, which     ,
     &                nev , np    , workl(irz),
     &                workl(bounds))
!
         if (msglvl .gt. 2) then
            call pzvout(comm,logfil, ncv, workl(irz), ndigit,
     &      '_neupd: Ritz values after calling _NGETS.')
            call pzvout(comm,logfil, ncv, workl(bounds), ndigit,
     &      '_neupd: Ritz value indices after calling _NGETS.')
         end if
!
!        %-----------------------------------------------------%
!        | Record indices of the converged wanted Ritz values  |
!        | Mark the select array for possible reordering       |
!        %-----------------------------------------------------%
!
         numcnv = 0
         do 11 j = 1,ncv
            rtemp = max(eps23,
     &                 dlapy2 ( dble (workl(irz+ncv-j)),
     &                          dimag(workl(irz+ncv-j)) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          dlapy2( dble (workl(ibd+jj-1)),
     &          dimag(workl(ibd+jj-1)) )
     &          .le. tol*rtemp) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
!
!        %-----------------------------------------------------------%
!        | Check the count (numcnv) of converged Ritz values with    |
!        | the number (nconv) reported by dnaupd.  If these two      |
!        | are different then there has probably been an error       |
!        | caused by incorrect passing of the dnaupd data.           |
!        %-----------------------------------------------------------%
!
         if (msglvl .gt. 2) then
             call pivout(comm, logfil, 1, [numcnv], ndigit,
     &            '_neupd: Number of specified eigenvalues')
             call pivout(comm, logfil, 1, [nconv], ndigit,
     &            '_neupd: Number of "converged" eigenvalues')
         end if
!
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
!
!        %-------------------------------------------------------%
!        | Call LAPACK routine zlahqr to compute the Schur form  |
!        | of the upper Hessenberg matrix returned by PZNAUPD.    |
!        | Make a copy of the upper Hessenberg matrix.           |
!        | Initialize the Schur vector matrix Q to the identity. |
!        %-------------------------------------------------------%
!
         call zcopy(ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call zlaset('All', ncv, ncv, zero, one, workl(invsub), ldq)
         call zlahqr(.true.       , .true.       , ncv          ,
     &                1            , ncv          , workl(iuptri),
     &                ldh          , workl(iheig) , 1            ,
     &                ncv          , workl(invsub), ldq          ,
     &                ierr         )
         call zcopy(ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
!
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
!
         if (msglvl .gt. 1) then
            call pzvout(comm, logfil, ncv, workl(iheig), ndigit,
     &           '_neupd: Eigenvalues of H')
            call pzvout(comm, logfil, ncv, workl(ihbds), ndigit,
     &           '_neupd: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call pzmout(comm, logfil, ncv, ncv,
     &              workl(iuptri), ldh, ndigit,
     &              '_neupd: The upper triangular matrix ')
            end if
         end if
         if (reord) then
!
!           %-----------------------------------------------%
!           | Reorder the computed upper triangular matrix. |
!           %-----------------------------------------------%
!
            call ztrsen('None'       , 'V'          , select      ,
     &                   ncv          , workl(iuptri), ldh         ,
     &                   workl(invsub), ldq          , workl(iheig),
     &                   nconv        , conds        , sep         ,
     &                   workev, ncv, ierr)
!
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
!
            if (msglvl .gt. 2) then
                call pzvout (comm, logfil, ncv, workl(iheig), ndigit,
     &           '_neupd: Eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call pzmout (comm, logfil, ncv, ncv,
     &                  workl(iuptri), ldq, ndigit,
     &              '_neupd: Triangular matrix after re-ordering')
                end if
            end if
         end if
!
!        %---------------------------------------------%
!        | Copy the last row of the Schur basis matrix |
!        | to workl(ihbds).  This vector will be used  |
!        | to compute the Ritz estimates of converged  |
!        | Ritz values.                                |
!        %---------------------------------------------%
!
         call zcopy(ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
!
!        %--------------------------------------------%
!        | Place the computed eigenvalues of H into D |
!        | if a spectral transformation was not used. |
!        %--------------------------------------------%
!
         if (type .eq. 'REGULR') then
            call zcopy(nconv, workl(iheig), 1, d, 1)
         end if
!
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(invsub,ldq).                            |
!        %----------------------------------------------------------%
!
         call zgeqr2(ncv, nconv , workl(invsub),
     &               ldq, workev, workev(ncv+1),
     &               ierr)
!
!        %--------------------------------------------------------%
!        | * Postmultiply V by Q using zunm2r.                    |
!        | * Copy the first NCONV columns of VQ into Z.           |
!        | * Postmultiply Z by R.                                 |
!        | The N by NCONV matrix Z is now a matrix representation |
!        | of the approximate invariant subspace associated with  |
!        | the Ritz values in workl(iheig). The first NCONV       |
!        | columns of V are now approximate Schur vectors         |
!        | associated with the upper triangular matrix of order   |
!        | NCONV in workl(iuptri).                                |
!        %--------------------------------------------------------%
!
         call zunm2r('Right', 'Notranspose', n            ,
     &                ncv    , nconv        , workl(invsub),
     &                ldq    , workev       , v            ,
     &                ldv    , workd(n+1)   , ierr         )
         call zlacpy('All', n, nconv, v, ldv, z, ldz)
!
         do 20 j=1, nconv
!
!           %---------------------------------------------------%
!           | Perform both a column and row scaling if the      |
!           | diagonal element of workl(invsub,ldq) is negative |
!           | I'm lazy and don't take advantage of the upper    |
!           | triangular form of workl(iuptri,ldq).             |
!           | Note that since Q is orthogonal, R is a diagonal  |
!           | matrix consisting of plus or minus ones.          |
!           %---------------------------------------------------%
!
            if ( dble( workl(invsub+(j-1)*ldq+j-1) ) .lt.
     &                  dble(zero) ) then
               call zscal(nconv, -one, workl(iuptri+j-1), ldq)
               call zscal(nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
!
 20      continue
!
         if (howmny .eq. 'A') then
!
!           %--------------------------------------------%
!           | Compute the NCONV wanted eigenvectors of T |
!           | located in workl(iuptri,ldq).              |
!           %--------------------------------------------%
!
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
!
            call ztrevc('Right', 'Select'     , select       ,
     &                   ncv    , workl(iuptri), ldq          ,
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , rwork        , ierr         )
!
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
!
!           %------------------------------------------------%
!           | Scale the returning eigenvectors so that their |
!           | Euclidean norms are all one. LAPACK subroutine |
!           | ztrevc returns each eigenvector normalized so  |
!           | that the element of largest magnitude has      |
!           | magnitude 1.                                   |
!           %------------------------------------------------%
!
            do 40 j=1, nconv
                  rtemp = dznrm2(ncv, workl(invsub+(j-1)*ldq), 1)
                  rtemp = dble(one) / rtemp
                  call zdscal ( ncv, rtemp,
     &                 workl(invsub+(j-1)*ldq), 1 )
!
!                 %------------------------------------------%
!                 | Ritz estimates can be obtained by taking |
!                 | the inner product of the last row of the |
!                 | Schur basis of H with eigenvectors of T. |
!                 | Note that the eigenvector matrix of T is |
!                 | upper triangular, thus the length of the |
!                 | inner product can be set to j.           |
!                 %------------------------------------------%
!
                  workev(j) = zzdotc(j, workl(ihbds), 1,
     &                        workl(invsub+(j-1)*ldq), 1)
 40         continue
!
            if (msglvl .gt. 2) then
               call zcopy(nconv, workl(invsub+ncv-1), ldq,
     &                    workl(ihbds), 1)
               call pzvout(comm, logfil, nconv, workl(ihbds), ndigit,
     &              '_neupd: Last row of the eigenvector matrix for T')
               if (msglvl .gt. 3) then
                  call pzmout(comm, logfil, nconv, ncv,
     &                 workl(invsub), ldq, ndigit,
     &                 '_neupd: The eigenvector matrix for T')
               end if
            end if
!
!           %---------------------------------------%
!           | Copy Ritz estimates into workl(ihbds) |
!           %---------------------------------------%
!
            call zcopy(nconv, workev, 1, workl(ihbds), 1)
!
!           %----------------------------------------------%
!           | The eigenvector matrix Q of T is triangular. |
!           | Form Z*Q.                                    |
!           %----------------------------------------------%
!
            call ztrmm('Right'   , 'Upper'      , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
!
         end if
!
      else
!
!        %--------------------------------------------------%
!        | An approximate invariant subspace is not needed. |
!        | Place the Ritz values computed PZNAUPD into D.    |
!        %--------------------------------------------------%
!
         call zcopy(nconv, workl(ritz), 1, d, 1)
         call zcopy(nconv, workl(ritz), 1, workl(iheig), 1)
         call zcopy(nconv, workl(bounds), 1, workl(ihbds), 1)
!
      end if
!
!     %------------------------------------------------%
!     | Transform the Ritz values and possibly vectors |
!     | and corresponding error bounds of OP to those  |
!     | of A*x = lambda*B*x.                           |
!     %------------------------------------------------%
!
      if (type .eq. 'REGULR') then
!
         if (rvec)
     &      call zscal(ncv, rnorm, workl(ihbds), 1)
!
      else
!
!        %---------------------------------------%
!        |   A spectral transformation was used. |
!        | * Determine the Ritz estimates of the |
!        |   Ritz values in the original system. |
!        %---------------------------------------%
!
         if (rvec)
     &      call zscal(ncv, rnorm, workl(ihbds), 1)
!
         do 50 k=1, ncv
            temp = workl(iheig+k-1)
            workl(ihbds+k-1) = workl(ihbds+k-1) / temp / temp
  50     continue
!
      end if
!
!     %-----------------------------------------------------------%
!     | *  Transform the Ritz values back to the original system. |
!     |    For TYPE = 'SHIFTI' the transformation is              |
!     |             lambda = 1/theta + sigma                      |
!     | NOTES:                                                    |
!     | *The Ritz vectors are not affected by the transformation. |
!     %-----------------------------------------------------------%
!
      if (type .eq. 'SHIFTI') then
         do 60 k=1, nconv
            d(k) = one / workl(iheig+k-1) + sigma
  60     continue
      end if
!
      if (type .ne. 'REGULR' .and. msglvl .gt. 1) then
         call pzvout (comm, logfil, nconv, d, ndigit,
     &     '_neupd: Untransformed Ritz values.')
         call pzvout (comm, logfil, nconv, workl(ihbds), ndigit,
     &     '_neupd: Ritz estimates of the untransformed Ritz values.')
      else if ( msglvl .gt. 1) then
         call pzvout (comm, logfil, nconv, d, ndigit,
     &     '_neupd: Converged Ritz values.')
         call pzvout (comm, logfil, nconv, workl(ihbds), ndigit,
     &     '_neupd: Associated Ritz estimates.')
      end if
!
!     %-------------------------------------------------%
!     | Eigenvector Purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 3. See reference 3.                  |
!     %-------------------------------------------------%
!
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
!
!        %------------------------------------------------%
!        | Purify the computed Ritz vectors by adding a   |
!        | little bit of the residual vector:             |
!        |                      T                         |
!        |          resid(:)*( e    s ) / theta           |
!        |                      NCV                       |
!        | where H s = s theta.                           |
!        %------------------------------------------------%
!
         do 100 j=1, nconv
            if (workl(iheig+j-1) .ne. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1)
     &                   /  workl(iheig+j-1)
            endif
 100     continue
!        %---------------------------------------%
!        | Perform a rank one update to Z and    |
!        | purify all the Ritz vectors together. |
!        %---------------------------------------%
!
         call zgeru(n, nconv, one, resid, 1, workev, 1, z, ldz)
!
      end if
!
 9000 continue
!
      return
!
!     %----------------%
!     | End of pzneupd |
!     %----------------%
!
      end
