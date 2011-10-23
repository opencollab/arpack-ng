#
# makefile to run simple examples of the reverse communication
# protocol.
#
# Modify if ARPACK library was built somewhere else.
# We assume that the required BLAS and LAPACK routines
# are in the version of libarpack.a built. If not, you
# will need to modify the link step below and link to them.
#
#\SCCS Information: @(#) 
# FILE: makefile   SID: 2.1   DATE OF SID: 11/19/95   RELEASE: 2
#
include ../../ARmake.inc
#
# Issue "make sym" to make all 12 symmetric drivers
# Issue "make ssdrv" to make 6 single precision symmetric drivers
# Issue "make dsdrv" to make 6 double precision symmetric drivers
#
sym: ssdrv dsdrv 
#
# simple symmetric problem using single precision
#
ssdrv: ssdrv1 ssdrv2 ssdrv3 ssdrv4 ssdrv5 ssdrv6

ssdrv1: ssdrv1.o 
	$(FC) $(FFLAGS) ssdrv1.o $(ALIBS) -o ssdrv1
#
ssdrv2: ssdrv2.o 
	$(FC) $(FFLAGS) ssdrv2.o $(ALIBS) -o ssdrv2
#
ssdrv3: ssdrv3.o 
	$(FC) $(FFLAGS) ssdrv3.o $(ALIBS) -o ssdrv3
#
ssdrv4: ssdrv4.o 
	$(FC) $(FFLAGS) ssdrv4.o $(ALIBS) -o ssdrv4
#
ssdrv5: ssdrv5.o 
	$(FC) $(FFLAGS) ssdrv5.o $(ALIBS) -o ssdrv5
#
ssdrv6: ssdrv6.o
	$(FC) $(FFLAGS) ssdrv6.o $(ALIBS) -o ssdrv6
#
#-----------------------------------------------------------------------
# simple symmetric problem using double precision
#
dsdrv: dsdrv1 dsdrv2 dsdrv3 dsdrv4 dsdrv5 dsdrv6

dsdrv1: dsdrv1.o
	$(FC) $(FFLAGS) dsdrv1.o $(ALIBS) -o dsdrv1
#
dsdrv2: dsdrv2.o 
	$(FC) $(FFLAGS) dsdrv2.o $(ALIBS) -o dsdrv2
#
dsdrv3: dsdrv3.o
	$(FC) $(FFLAGS) dsdrv3.o $(ALIBS) -o dsdrv3
#
dsdrv4: dsdrv4.o
	$(FC) $(FFLAGS) dsdrv4.o $(ALIBS) -o dsdrv4
#
dsdrv5: dsdrv5.o
	$(FC) $(FFLAGS) dsdrv5.o $(ALIBS) -o dsdrv5
#
dsdrv6: dsdrv6.o
	$(FC) $(FFLAGS) dsdrv6.o $(ALIBS) -o dsdrv6 
