1. Purpose 
   -------
   This directory contains example drivers that call ARPACK subroutines
   [s,d]saupd.f and [s,d]seupd.f to solve SYMMETRIC eigenvalue problems 
   using regular, inverse, shift-invert or other special modes (such as 
   Cayley, Bucking etc.)  

   These drivers illustrate how to set various ARPACK parameters to solve 
   different problems in different modes.  They provide a guideline on how 
   to use ARPACK's reverse communication interface.  The user may modify
   any one of these drivers, and provide his/her own matrix vector
   multiplication routine to solve the problem of his/her own interest.


2. Naming convention
   -----------------
  
   The name for each driver has the form 'XsdrvN.f', where
   X - is 's' (single precision)
       or 'd' (double precision)

   N - is a number between 1 and 6.  If 
       N = 1, the driver solves a STANDARD eigenvalue problem in
              REGULAR mode
       N = 2, the driver solves a STANDARD eigenvalue problem in
              SHIFT-INVERT mode. 
       N = 3, the driver solves a GENERALIZED eigenvalue problem in
              INVERSE mode
       N = 4, the driver solves a GENERALIZED eigenvalue problem in
              SHIFT-INVERT mode.

   These are 4 commonly used drivers.  For shift-invert (N=2,4) mode
   the user needs to supply a linear system solver to perform
   y=inv[A-sigma*B]*x.

   When N > 4, a special mode is used.  If  
       N = 5, the driver solves a GENERALIZED eigenvalue problem in
              BUCKLING mode.
       N = 6, the driver solves a GENERALIZED eigenvalue problem in
              CAYLEY mode.

   These two drivers require the user to provide linear system
   solvers also.  For more information on Cayley and Buckling mode,
   see the following reference:

   R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
   Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
   SIAM J. Matr. Anal. Apps.,  January (1993).


3. Usage
   -----
   To run these drivers, you may use the makefile in this
   directory and issue, for example, "make ssdrv1".  Then
   execute using "ssdrv1".
