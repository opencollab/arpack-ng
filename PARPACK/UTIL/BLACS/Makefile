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
#  The name of the library is defined by $(ARPACKLIB) in
#  ../ARmake.inc and is created at the next higher directory level.
#

OBJS  = pivout.o

SOBJ  = psvout.o  psmout.o

DOBJ  = pdvout.o  pdmout.o

COBJ  = pcvout.o  pcmout.o 

ZOBJ  = pzvout.o  pzmout.o

.SUFFIXES:      .o .F .f

.f.o:
	$(FC) $(FFLAGS) -c $<
#
#  make the library containing both single and double precision
#
all: single complex double complex16
 
single: $(SOBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(SOBJ) $(OBJS)
	$(RANLIB) $(PARPACKLIB)
 
double: $(DOBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(DOBJ) $(OBJS)
	$(RANLIB) $(PARPACKLIB)
 
complex: $(COBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(COBJ) $(OBJS)
	$(RANLIB) $(PARPACKLIB)
 
complex16: $(ZOBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(PARPACKLIB) $(ZOBJ) $(OBJS)
	$(RANLIB) $(PARPACKLIB)

#
#  clean	- remove all object files
#
clean:
	rm -f *.o a.out core

