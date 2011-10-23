1. Purpose 
   -------
   This directory contains example drivers that call subroutines
   [s,d]sband.f, [s,d]nband.f or [c,z]nband.f to solve various eigenvalue 
   problems in which matrices have BANDED structure. These drivers illustrate 
   how to construct LAPACK-style banded matrices, and how to set various 
   ARPACK parameters to solve different problems in different modes. The 
   user may modify any one of these drivers, and construct his/her own matrix 
   to solve the problem of his/her own interest.

2. Naming convention
   -----------------
   The name for each banded solver has the form 'XYband.f', where
   X - is 's' (single precision)
       or 'd' (double precision)
       or 'c' (single precision complex)
       or 'z' (double precision complex)

   Y - is 's' (symmetric)
       or 'n' (nonsymmetric)

   The name for each driver has the form 'XYbdrN.f', where
   X - is 's' (single precision)
       or 'd' (double precision)
       or 'c' (single precision complex)
       or 'z' (double precision complex)

   Y - is 's' (symmetric)
       or 'n' (nonsymmetric)

   N - is a number between 1 and 6.  If
       N = 1, the driver solves a STANDARD eigenvalue problem in
              REGULAR mode.
       N = 2, the driver solves a STANDARD eigenvalue problem in
              SHIFT-INVERT mode.
       N = 3, the driver solves a GENERALIZED eigenvalue problem in
              INVERSE mode.
       N = 4, the driver solves a GENERALIZED eigenvalue problem in
              SHIFT-INVERT mode (using mode 3 of __aupd.)
   These are four commonly used drivers.

   When N > 4 (only for real matrices), a special mode is used.
   For symmetric problem, if
       N = 5, the driver solves a GENERALIZED eigenvalue problem in
              BUCKLING mode.
       N = 6, the driver solves a GENERALIZED eigenvalue problem in
              CAYLEY mode.

   For nonsymmetric problem, if
       N = 5, the driver solves a STANDARD eigenvalue problem in
              SHIFT-INVERT mode using mode 4 of [d,s]naupd.f.
       N = 6. the driver solves a GENERALIZED eigenvalue problem
              SHIFT-INVERT mode using mode 4 of [d,s]naupd.f.
       Note: the imaginary part of the shift MUST be nonzero when
       these two drivers are used.

3. Usage
   -----
   To run these drivers, you may use the makefile in this
   directory and issue, for example, "make snbdr1".  Then
   execute using "snbdr1".
 
