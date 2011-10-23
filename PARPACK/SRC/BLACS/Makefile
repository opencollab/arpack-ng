############################################################################
#
#  Program:         PARPACK
#
#  Module:          Makefile
#
#  Purpose:         Sources Makefile
#
#  Creation date:   February 22, 1996
#
#  Modified:
#
#  Send bug reports, comments or suggestions to arpack.caam.rice.edu
#
############################################################################

include ../../../ARmake.inc

############################################################################
#  To create or add to the library, enter make followed by one or
#  more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is defined by $(PARPACKLIB) in
#  $(home)/ARmake.inc and is created in the $(home) directory.

SOBJ  = psgetv0.o  \
	psnaitr.o psnapps.o psnaup2.o psnaupd.o psneigh.o psngets.o \
	pssaitr.o pssapps.o pssaup2.o pssaupd.o psseigt.o pssgets.o \
	psneupd.o psseupd.o pslarnv.o pslamch.o psnorm2.o

DOBJ  = pdgetv0.o \
	pdnaitr.o pdnapps.o pdnaup2.o pdnaupd.o pdneigh.o pdngets.o \
	pdsaitr.o pdsapps.o pdsaup2.o pdsaupd.o pdseigt.o pdsgets.o \
	pdneupd.o pdseupd.o pdlarnv.o pdlamch.o pdnorm2.o

COBJ  = pcnaitr.o pcnapps.o pcnaup2.o pcnaupd.o pcneigh.o \
		pcneupd.o pcngets.o pcgetv0.o pscnorm2.o pclarnv.o

ZOBJ  = pznaitr.o pznapps.o pznaup2.o pznaupd.o pzneigh.o \
		pzneupd.o pzngets.o pzgetv0.o pdznorm2.o pzlarnv.o

.f.o:
	$(FC) $(FFLAGS) -c $<

all: single complex double complex16

single: $(SOBJ)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(SOBJ)
	$(RANLIB) $(PARPACKLIB)

double: $(DOBJ)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(DOBJ)
	$(RANLIB) $(PARPACKLIB)

complex: $(COBJ)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(COBJ)
	$(RANLIB) $(PARPACKLIB)

complex16: $(ZOBJ)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(ZOBJ)
	$(RANLIB) $(PARPACKLIB)

#
#  clean	- remove all object files
#
clean:
	rm -f *.o a.out core

