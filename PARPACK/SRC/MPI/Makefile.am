F77 = $(MPIF77)
noinst_LTLIBRARIES = libparpack_noopt.la libparpacksrcmpi.la

FFLAGS_SAV = @FFLAGS@
FFLAGS = 
libparpack_noopt_la_SOURCES = pslamch.f pdlamch.f
libparpack_noopt_la_FFLAGS = -O0

libparpacksrcmpi_la_SOURCES = \
	psgetv0.f  \
	psnaitr.f psnapps.f psnaup2.f psnaupd.f psneigh.f psngets.f \
	pssaitr.f pssapps.f pssaup2.f pssaupd.f psseigt.f pssgets.f \
	psneupd.f psseupd.f pslarnv.f psnorm2.f pdgetv0.f \
	pdnaitr.f pdnapps.f pdnaup2.f pdnaupd.f pdneigh.f pdngets.f \
	pdsaitr.f pdsapps.f pdsaup2.f pdsaupd.f pdseigt.f pdsgets.f \
	pdneupd.f pdseupd.f pdlarnv.f pdnorm2.f \
	pcnaitr.f pcnapps.f pcnaup2.f pcnaupd.f pcneigh.f \
	pcneupd.f pcngets.f pcgetv0.f pscnorm2.f pclarnv.f \
	pznaitr.f pznapps.f pznaup2.f pznaupd.f pzneigh.f \
	pzneupd.f pzngets.f pzgetv0.f pdznorm2.f pzlarnv.f

libparpacksrcmpi_la_FFLAGS = $(FFLAGS_SAV)
libparpacksrcmpi_la_LIBADD = libparpack_noopt.la
