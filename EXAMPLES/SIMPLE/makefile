#
# makefile to run simple examples of the reverse communication
# protocol.
#
# Modify if ARPACK library was built somewhere else.
#
#\SCCS Information: @(#) 
# FILE: makefile   SID: 2.2   DATE OF SID: 9/24/96   RELEASE: 2
#
include ../../ARmake.inc
#
# simple drivers
#
simple: sssimp dssimp snsimp dnsimp cnsimp znsimp
#
sssimp: sssimp.o
	$(FC) $(FFLAGS) sssimp.o $(ALIBS) -o sssimp
dssimp: dssimp.o
	$(FC) $(FFLAGS) dssimp.o $(ALIBS) -o dssimp
snsimp: snsimp.o
	$(FC) $(FFLAGS) snsimp.o $(ALIBS) -o snsimp
dnsimp: dnsimp.o
	$(FC) $(FFLAGS) dnsimp.o $(ALIBS) -o dnsimp
cnsimp: cnsimp.o
	$(FC) $(FFLAGS) cnsimp.o $(ALIBS) -o cnsimp 
znsimp: znsimp.o
	$(FC) $(FFLAGS) znsimp.o $(ALIBS) -o znsimp  
