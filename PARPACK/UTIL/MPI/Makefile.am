F77 = $(MPIF77)

SRCS = pivout.f

SSRC = psmout.f psvout.f
DSRC = pdmout.F90 pdvout.f
CSRC = pcmout.F90 pcvout.F90
ZSRC = pzmout.f pzvout.f

noinst_LTLIBRARIES = libparpackutil.la
libparpackutil_la_SOURCES = $(SRCS) $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
