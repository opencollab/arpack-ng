#
# makefile to run simple examples of the reverse communication
# protocol.
#
# Modify if ARPACK library was built somewhere else.
#
# If ALIBS (defined in ../../ARmake.inc) contains the BLAS and LAPACK
# libraries installed on your system, you DO NOT need to change this
# makefile. OTHERWISE, you may need to modify the Makefile in the top level
# ARPACK directory tree to include sdrv and/or ddrv in the
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
# Issue "make nonsym" to make all 12 nonsymmetric drivers.
# Issue "make sndrv" to make 6 single precision nonsymmetric drivers.
# Issue "make dndrv" to make 6 double precision nonsymmetric drivers.
#
nonsym: sndrv dndrv 
#
#-----------------------------------------------------------------------
# Simple nonsymmetric problem using single precision
#
sndrv: sndrv1 sndrv2 sndrv3 sndrv4 sndrv5 sndrv6

sndrv1: sndrv1.o  
	$(FC) $(FFLAGS) sndrv1.o $(ALIBS)  -o sndrv1
#
sndrv2: sndrv2.o 
	$(FC) $(FFLAGS) sndrv2.o $(ALIBS) -o sndrv2
#
sndrv3: sndrv3.o 
	$(FC) $(FFLAGS) sndrv3.o $(ALIBS) -o sndrv3
#
sndrv4: sndrv4.o 
	$(FC) $(FFLAGS) sndrv4.o $(ALIBS) -o sndrv4
#
sndrv5: sndrv5.o 
	$(FC) $(FFLAGS) sndrv5.o $(ALIBS) -o sndrv5
#
sndrv6: sndrv6.o 
	$(FC) $(FFLAGS) sndrv6.o $(ALIBS) -o sndrv6
#
#------------------------------------------------------
# Simple nonsymmetric problem using double precision
#
dndrv: dndrv1 dndrv2 dndrv3 dndrv4 dndrv5 dndrv6

dndrv1: dndrv1.o 
	$(FC) $(FFLAGS) dndrv1.o $(ALIBS) -o dndrv1 
#
dndrv2: dndrv2.o 
	$(FC) $(FFLAGS) dndrv2.o $(ALIBS) -o dndrv2
#
dndrv3: dndrv3.o 
	$(FC) $(FFLAGS) dndrv3.o $(ALIBS) -o dndrv3
#
dndrv4: dndrv4.o 
	$(FC) $(FFLAGS) dndrv4.o $(ALIBS) -o dndrv4
#
dndrv5: dndrv5.o 
	$(FC) $(FFLAGS) dndrv5.o $(ALIBS) -o dndrv5
#
dndrv6: dndrv6.o 
	$(FC) $(FFLAGS) dndrv6.o $(ALIBS) -o dndrv6
