1. Purpose 
   -------
   This directory contains example drivers that call ARPACK subroutines
   [c,z]naupd.f and [c,z]neupd.f to solve COMPLEX eigenvalue problems 
   using regular, inverse or shift-invert modes. These drivers illustrate 
   how to set various ARPACK parameters to solve different problems in 
   different modes.  They provide a guideline on how to use ARPACK's reverse 
   communication interface.  The user may modify any one of these drivers, 
   and provide his/her own matrix vector multiplication routine to solve 
   the problem of his/her own interest.

2. Naming convention
   -----------------

   The name for each driver has the form 'XndrvN.f', where
   X - is 'c' (single precision complex)
       or 'z' (double precision complex)

   N - is a number between 1 and 4.  If
       N = 1, the driver solves a STANDARD eigenvalue problem in
              REGULAR mode
       N = 2, the driver solves a STANDARD eigenvalue problem in
              SHIFT-INVERT mode.
       N = 3, the driver solves a GENERALIZED eigenvalue problem in
              INVERSE mode
       N = 4, the driver solves a GENERALIZED eigenvalue problem in
              SHIFT-INVERT mode.

   For shift-invert modes (N=2,4), the user needs to provide a complex
   linear system solver to perform y=inv[A-sigma*B]*x.

3. Usage
   -----
   To run these drivers, you may use the makefile in this
   directory and issue, for example, "make cndrv1".  Then
   execute using "cndrv1".
