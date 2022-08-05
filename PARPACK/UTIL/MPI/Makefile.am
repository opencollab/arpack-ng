F77 = $(MPIF77)

SRCS = pivout.F90

SSRC = psmout.F90 psvout.F90
DSRC = pdmout.F90 pdvout.F90
CSRC = pcmout.F90 pcvout.F90
ZSRC = pzmout.f pzvout.f

noinst_LTLIBRARIES = libparpackutil.la
libparpackutil_la_SOURCES = $(SRCS) $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
