!\BeginDoc
!
!\Name: pzlarnv
!
! Message Passing Layer: MPI
!
!\Description:
!
!  Parallel Version of ARPACK utility routine zlarnv
!
!  PZLARNV returns a vector of n (nloc) random Complex*16 numbers from a uniform or
!  normal distribution. It is assumed that X is distributed across a 1-D array
!  of processors ( nprocs < 1000 )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated.
!
!  X       (output) Complex*16 array, dimension (N)
!          The generated random numbers.
!
!\Author: Kristi Maschhoff
!
!\Details
!
!  Simple parallel version of LAPACK auxiliary routine zlarnv
!  for X distributed across a 1-D array of processors.
!  This routine calls the auxiliary routine CLARNV to generate random
!  Complex*16 numbers from a uniform or normal distribution. Output is consistent
!  with serial version.
!
!\SCCS Information:
! FILE: larnv.F   SID: 1.3   DATE OF SID: 04/17/99
!
!-----------------------------------------------------------------------
!
      subroutine pzlarnv( comm, idist, iseed, n, x )
!
      integer   comm
!     ..
!     .. Scalar Arguments ..
      integer			idist, n
!     ..
!     .. Array Arguments ..
      integer			iseed( 4 )
      Complex*16
     &                          x( * )
!     ..
!     .. External Subroutines ..
      external			zlarnv
!     ..
!     .. Executable Statements ..
!
      call zlarnv ( idist, iseed, n, x )
!
      return
      end
