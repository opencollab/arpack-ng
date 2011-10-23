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

ssvd: ssvd.o
	$(FC) $(FFLAGS) ssvd.o $(ALIBS) -o ssvd
dsvd: dsvd.o
	$(FC) $(FFLAGS) dsvd.o $(ALIBS) -o dsvd
