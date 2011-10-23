1. Purpose 
   -------
   This directory contains example drivers that call ARPACK subroutines
   [s,d]naupd.f and [s,d]neupd.f to solve real NONSYMMETRIC eigenvalue 
   problems using regular, inverse or shift-invert modes.  These drivers 
   illustrate how to set various ARPACK parameters to solve different 
   problems in different modes.  They provide a guideline on how to use 
   ARPACK's reverse communication interface.  The user may modify any one 
   of these drivers, and provide his/her own matrix vector multiplication 
   routine to solve the problem of his/her own interest.

2. Naming convention
   -----------------
   The name for each driver has the form 'XndrvN.f', where
   X - is 's' (single precision)
       or 'd' (double precision)

   N - is a number between 1 and 6.  If
       N = 1, the driver solves a STANDARD eigenvalue problem in
              REGULAR mode.
       N = 2, the driver solves a STANDARD eigenvalue problem in
              SHIFT-INVERT mode with a REAL shift.
       N = 3, the driver solves a GENERALIZED eigenvalue problem in
              INVERSE mode.
       N = 4, the driver solves a GENERALIZED eigenvalue problem in
              SHIFT-INVERT mode with a REAL shift.

   These are 4 commonly used drivers.  For shift-invert (N=2,4) mode
   the user needs to supply a linear system solver to perform
   y=inv[A-sigma*B]*x.

   If N > 4, shift-invert is used with a complex shift whose imaginary 
   part is nonzero. If
       N = 5, the driver solves a GENERALIZED eigenvalue problem 
              using mode 3 of [s,d]naupd. 
       N = 6. the driver solves a GENERALIZED eigenvalue problem 
              using mode 4 of [s,d]naupd.

   These two drivers require the user to provide COMPLEX arithmetic
   linear system solver.  For more information on the use of
   complex shift, see the following reference:

   B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
   Real Matrices", Linear Algebra and its Applications, vol 88/89,
   pp 575-595, (1987).

3. Usage
   -----
   To run these drivers, you may use the makefile in this
   directory and issue, for example, "make sndrv1".  Then
   execute using "sndrv1".
