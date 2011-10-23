#
# makefile to run simple examples of the reverse communication
# protocol.
#
# Modify if ARPACK library was built somewhere else.
# We assume that the required BLAS and LAPACK routines
# are in the version of libarpack.a built. If not, you
# will need to modify the link step below and link to them.
#
#
include ../../../ARmake.inc

# Issue "make nonsym" to make all nonsymmetric drivers.
# Issue "make psndrv" to make single precision nonsymmetric drivers.
# Issue "make pdndrv" to make double precision nonsymmetric drivers.
#
nonsym: psndrv pdndrv 
#
# Issue "make sym" to make all symmetric drivers
# Issue "make pssdrv" to make single precision symmetric drivers
# Issue "make pdsdrv" to make double precision symmetric drivers
#
sym: pssdrv pdsdrv
#
# Issue "complex" to make all complex drivers.
# Issue "pcndrv" to make single precision complex drivers.
# Issue "pzndrv" to make double precision complex drivers.
#
complex: pcndrv pzndrv
#
#-----------------------------------------------------------------------
# Simple nonsymmetric problem using single precision
#
psndrv: psndrv1 psndrv3 

psndrv1: psndrv1.o  
	$(PFC) $(PFFLAGS) psndrv1.o $(PLIBS) -o psndrv1_$(PLAT)
#
psndrv3: psndrv3.o 
	$(PFC) $(PFFLAGS) psndrv3.o $(PLIBS) -o psndrv3_$(PLAT)
#
#------------------------------------------------------
# Simple nonsymmetric problem using double precision
#
pdndrv: pdndrv1 pdndrv3 

pdndrv1: pdndrv1.o 
	$(PFC) $(PFFLAGS) pdndrv1.o $(PLIBS) -o pdndrv1_$(PLAT)
#
pdndrv3: pdndrv3.o 
	$(PFC) $(PFFLAGS) pdndrv3.o $(PLIBS) -o pdndrv3_$(PLAT)
#
#-----------------------------------------------------------------------
# Simple symmetric problem using single precision
#
pssdrv: pssdrv1 

pssdrv1: pssdrv1.o
	$(PFC) $(PFFLAGS) pssdrv1.o $(PLIBS) -o pssdrv1_$(PLAT)
#
#------------------------------------------------------
# Simple symmetric problem using double precision
#
pdsdrv: pdsdrv1 

pdsdrv1: pdsdrv1.o
	$(PFC) $(PFFLAGS) pdsdrv1.o $(PLIBS) -o pdsdrv1_$(PLAT)

#-----------------------------------------------------------------------
# Complex problem using single complex
#
pcndrv: pcndrv1

pcndrv1: pcndrv1.o
	$(PFC) $(PFFLAGS) pcndrv1.o $(PLIBS) -o pcndrv1_$(PLAT)
#
#----------------------------------------------------------------------
# Complex problem using double complex
#
pzndrv: pzndrv1

pzndrv1: pzndrv1.o
	$(PFC) $(PFFLAGS) pzndrv1.o $(PLIBS) -o pzndrv1_$(PLAT)
