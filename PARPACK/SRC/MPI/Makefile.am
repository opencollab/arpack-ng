F77 = $(MPIF77)
FFLAGS_SAV = @FFLAGS@
FFLAGS =

SSRC = psnaitr.f psnapps.f psnaup2.f psnaupd.f psneigh.f psneupd.f psngets.f \
       pssaitr.f pssapps.f pssaup2.f pssaupd.f psseigt.f psseupd.f pssgets.f \
       psgetv0.f pslarnv.f psnorm2.f

DSRC = pdnaitr.f pdnapps.f pdnaup2.f pdnaupd.f pdneigh.f pdneupd.f pdngets.f \
       pdsaitr.f pdsapps.f pdsaup2.f pdsaupd.f pdseigt.f pdseupd.f pdsgets.f \
       pdgetv0.f pdlarnv.f pdnorm2.f

CSRC = pcnaitr.f pcnapps.f pcnaup2.f pcnaupd.f pcneigh.f pcneupd.f pcngets.f \
       pcgetv0.f pclarnv.f pscnorm2.f

ZSRC = pznaitr.f pznapps.f pznaup2.f pznaupd.f pzneigh.f pzneupd.f pzngets.f \
       pzgetv0.f pzlarnv.f pdznorm2.f

EXTRA_DIST = debug.h stat.h

noinst_LTLIBRARIES = libparpack_noopt.la libparpacksrc.la

libparpack_noopt_la_SOURCES = pslamch.f pdlamch10.f
libparpack_noopt_la_FFLAGS = -O0

libparpacksrc_la_SOURCES = $(SSRC) $(DSRC) $(CSRC) $(ZSRC)
libparpacksrc_la_FFLAGS = $(FFLAGS_SAV)
libparpacksrc_la_LIBADD = libparpack_noopt.la
