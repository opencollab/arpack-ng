F77 = $(MPIF77)

SRCS = pivout.f

SSRC = psmout.f psvout.f
DSRC = pdmout.f pdvout.f
CSRC = pcmout.f pcvout.f
ZSRC = pzmout.f pzvout.f

noinst_LTLIBRARIES = libparpackutil.la
libparpackutil_la_SOURCES = $(SRCS) $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
