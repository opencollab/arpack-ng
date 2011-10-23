#
# makefile to run simple examples of the reverse communication
# protocol.
#
# Modify if ALIBS library was built somewhere else.
#
# If ALIBS (defined in ../../ARmake.inc) contains the BLAS and LAPACK
# libraries installed on your system, you DO NOT need to change this
# makefile. OTHERWISE, you may need to modify the Makefile in the top level
# ARPACK directory tree to include cdrv and/or zdrv in the
# definition of the PRECISION variable, and issue 'make lib' from there.
# For example:
#
#         PRECISION = single double sdrv ddrv
#
#\SCCS Information: @(#) 
# FILE: makefile   SID: 2.2   DATE OF SID: 9/24/96   RELEASE: 2 
#
include ../../ARmake.inc
#
# Issue "complex" to make all 8 complex drivers.
# Issue "cndrv" to make 4 single precision complex drivers.
# Issue "zndrv" to make 4 double precision complex drivers.
#

complex: cndrv zndrv 

#-----------------------------------------------------------------------
# Complex problem using single complex
#
cndrv: cndrv1 cndrv2 cndrv3 cndrv4

cndrv1: cndrv1.o 
	$(FC) $(FFLAGS) cndrv1.o $(ALIBS) -o cndrv1
#
cndrv2: cndrv2.o 
	$(FC) $(FFLAGS) cndrv2.o $(ALIBS) -o cndrv2
#
cndrv3: cndrv3.o 
	$(FC) $(FFLAGS) cndrv3.o $(ALIBS) -o cndrv3
#
cndrv4: cndrv4.o 
	$(FC) $(FFLAGS) cndrv4.o $(ALIBS) -o cndrv4
#
#----------------------------------------------------------------------
# Complex problem using double complex
#
zndrv: zndrv1 zndrv2 zndrv3 zndrv4

zndrv1: zndrv1.o 
	$(FC) $(FFLAGS) zndrv1.o $(ALIBS) -o zndrv1
#
zndrv2: zndrv2.o
	$(FC) $(FFLAGS) zndrv2.o $(ALIBS) -o zndrv2
#
zndrv3: zndrv3.o 
	$(FC) $(FFLAGS) zndrv3.o $(ALIBS) -o zndrv3
#
zndrv4: zndrv4.o 
	$(FC) $(FFLAGS) zndrv4.o $(ALIBS) -o zndrv4
