#
# makefile to run simple examples of the reverse communication
# protocol.
#
# Modify if ARPACK library was built somewhere else.
#
# If ALIBS (defined in ../../ARmake.inc) contains the BLAS and LAPACK
# libraries installed on your system, you DO NOT need to change this
# makefile. OTHERWISE, you may need to modify the Makefile in the top level
# ARPACK directory tree to include sdrv, ddrv, cdrv, or zdrv in the
# definition of the PRECISION variable, and issue 'make lib' there.
# For example:
#
#         PRECISION = single double sdrv ddrv
#
#\SCCS Information: @(#) 
# FILE: makefile   SID: 2.2   DATE OF SID: 9/24/96   RELEASE: 2 
#
include ../../ARmake.inc

# Issue "make band" to make all 32 band drivers.
# Issue "make ssbdr" to make 6 single precision symmetric band drivers.
# Issue "make dsbdr" to make 6 double precision symmetric band drivers.
# Issue "make snbdr" to make 6 single precision nonsymmetric band drivers.
# Issue "make dnbdr" to make 6 double precision nonsymmetric band drivers.
# Issue "make cnbdr" to make 4 single precision complex band drivers.
# Issue "make znbdr" to make 4 double precision complex band drivers. 
#

band: ssbdr dsbdr snbdr dnbdr cnbdr znbdr

#
#---------------------------------------------------------------------
# Simple symetric problem using BAND solver (single precision)
#
ssbdr: ssbdr1 ssbdr2 ssbdr3 ssbdr4 ssbdr5 ssbdr6

ssbdr1: ssbdr1.o ssband.o
	$(FC) $(FFLAGS) ssbdr1.o ssband.o $(ALIBS) -o ssbdr1
#
ssbdr2: ssbdr2.o ssband.o
	$(FC) $(FFLAGS) ssbdr2.o ssband.o $(ALIBS) -o ssbdr2
#
ssbdr3: ssbdr3.o ssband.o
	$(FC) $(FFLAGS) ssbdr3.o ssband.o $(ALIBS) -o ssbdr3
#
ssbdr4: ssbdr4.o ssband.o
	$(FC) $(FFLAGS) ssbdr4.o ssband.o $(ALIBS) -o ssbdr4
#
ssbdr5: ssbdr5.o ssband.o
	$(FC) $(FFLAGS) ssbdr5.o ssband.o $(ALIBS) -o ssbdr5
#
ssbdr6: ssbdr6.o ssband.o
	$(FC) $(FFLAGS) ssbdr6.o ssband.o $(ALIBS) -o ssbdr6
#
#-------------------------------------------------------------
#
# Simple symmetric problem using BAND solver (double precision)
#
dsbdr: dsbdr1 dsbdr2 dsbdr3 dsbdr4 dsbdr5 dsbdr6

dsbdr1: dsbdr1.o dsband.o
	$(FC) $(FFLAGS) dsbdr1.o dsband.o $(ALIBS) -o dsbdr1
#
dsbdr2: dsbdr2.o dsband.o 
	$(FC) $(FFLAGS) dsbdr2.o dsband.o $(ALIBS) -o dsbdr2
#
dsbdr3: dsbdr3.o dsband.o
	$(FC) $(FFLAGS) dsbdr3.o dsband.o $(ALIBS) -o dsbdr3
#
dsbdr4: dsbdr4.o dsband.o
	$(FC) $(FFLAGS) dsbdr4.o dsband.o $(ALIBS) -o dsbdr4
#
dsbdr5: dsbdr5.o dsband.o
	$(FC) $(FFLAGS) dsbdr5.o dsband.o $(ALIBS) -o dsbdr5
#
dsbdr6: dsbdr6.o dsband.o
	$(FC) $(FFLAGS) dsbdr6.o dsband.o $(ALIBS) -o dsbdr6
#
#----------------------------------------------------------------------
# Simple nonsymmetric problem using BAND solver (single precision)
#
snbdr: snbdr1 snbdr2 snbdr3 snbdr4 snbdr5 snbdr6

snbdr1: snbdr1.o snband.o 
	$(FC) $(FFLAGS) snbdr1.o snband.o $(ALIBS) -o snbdr1
#
snbdr2: snbdr2.o snband.o 
	$(FC) $(FFLAGS) snbdr2.o snband.o $(ALIBS) -o snbdr2
#
snbdr3: snbdr3.o snband.o 
	$(FC) $(FFLAGS) snbdr3.o snband.o $(ALIBS) -o snbdr3
#
snbdr4: snbdr4.o snband.o 
	$(FC) $(FFLAGS) snbdr4.o snband.o $(ALIBS) -o snbdr4
#
snbdr5: snbdr5.o snband.o 
	$(FC) $(FFLAGS) snbdr5.o snband.o $(ALIBS) -o snbdr5
#
snbdr6: snbdr6.o snband.o 
	$(FC) $(FFLAGS) snbdr6.o snband.o $(ALIBS) -o snbdr6
#
#-------------------------------------------------------------------------
# Simple nonsymmetric problem using BAND solver (double precision)
#
dnbdr: dnbdr1 dnbdr2 dnbdr3 dnbdr4 dnbdr5 dnbdr6

dnbdr1: dnbdr1.o dnband.o 
	$(FC) $(FFLAGS) dnbdr1.o dnband.o $(ALIBS) -o dnbdr1
#
dnbdr2: dnbdr2.o dnband.o 
	$(FC) $(FFLAGS) dnbdr2.o dnband.o $(ALIBS) -o dnbdr2
#
dnbdr3: dnbdr3.o dnband.o 
	$(FC) $(FFLAGS) dnbdr3.o dnband.o $(ALIBS) -o dnbdr3
#
dnbdr4: dnbdr4.o dnband.o 
	$(FC) $(FFLAGS) dnbdr4.o dnband.o $(ALIBS) -o dnbdr4
#
dnbdr5: dnbdr5.o dnband.o 
	$(FC) $(FFLAGS) dnbdr5.o dnband.o $(ALIBS) -o dnbdr5
#
dnbdr6: dnbdr6.o dnband.o 
	$(FC) $(FFLAGS) dnbdr6.o dnband.o $(ALIBS) -o dnbdr6
#
#-------------------------------------------------------------------------
# Complex problem using BAND solver (single precision complex)
#
cnbdr:  cnbdr1 cnbdr2 cnbdr3 cnbdr4

cnbdr1: cnbdr1.o cnband.o
	$(FC) $(FFLAGS) cnbdr1.o cnband.o $(ALIBS) -o cnbdr1
#
cnbdr2: cnbdr2.o cnband.o
	$(FC) $(FFLAGS) cnbdr2.o cnband.o $(ALIBS) -o cnbdr2
#
cnbdr3: cnbdr3.o cnband.o
	$(FC) $(FFLAGS) cnbdr3.o cnband.o $(ALIBS) -o cnbdr3
#
cnbdr4: cnbdr4.o cnband.o
	$(FC) $(FFLAGS) cnbdr4.o cnband.o $(ALIBS) -o cnbdr4
#
#-------------------------------------------------------------------------
# Complex problem using BAND solver (double precision complex)
#
znbdr: znbdr1 znbdr2 znbdr3 znbdr4

znbdr1: znbdr1.o znband.o
	$(FC) $(FFLAGS) znbdr1.o znband.o $(ALIBS) -o znbdr1
#
znbdr2: znbdr2.o znband.o
	$(FC) $(FFLAGS) znbdr2.o znband.o $(ALIBS) -o znbdr2
#
znbdr3: znbdr3.o znband.o
	$(FC) $(FFLAGS) znbdr3.o znband.o $(ALIBS) -o znbdr3
#
znbdr4: znbdr4.o znband.o
	$(FC) $(FFLAGS) znbdr4.o znband.o $(ALIBS) -o znbdr4
