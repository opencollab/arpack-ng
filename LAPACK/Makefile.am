noinst_LTLIBRARIES = liblapack_noopt.la liblapack_arpack.la

# liblapack_noopt_la_FFLAGS will not overwrite user flags
FFLAGS_SAV = @FFLAGS@
FFLAGS =
liblapack_noopt_la_SOURCES = \
	slamch.f slaruv.f dlamch.f dlaruv.f
liblapack_noopt_la_FFLAGS = -O0

liblapack_arpack_la_SOURCES = \
	sgeqr2.f slabad.f slacon.f slacpy.f sladiv.f slae2.f slaev2.f \
        slaexc.f slagtm.f slahqr.f slaln2.f slange.f slanhs.f \
        slanst.f slanv2.f slaptm.f slapy2.f slapy3.f slaran.f slarf.f \
        slarfg.f slarfx.f slarnd.f slarnv.f slartg.f slascl.f \
        slaset.f slasr.f  slasrt.f slassq.f slasy2.f sorm2r.f ssteqr.f \
        strevc.f strexc.f strsen.f strsyl.f \
	dgeqr2.f dlabad.f dlacon.f dlacpy.f dladiv.f dlae2.f dlaev2.f \
        dlaexc.f dlagtm.f dlahqr.f dlaln2.f dlange.f dlanhs.f \
        dlanst.f dlanv2.f dlaptm.f dlapy2.f dlapy3.f dlaran.f dlarf.f \
        dlarfg.f dlarfx.f dlarnd.f dlarnv.f dlartg.f dlascl.f \
        dlaset.f dlasr.f  dlasrt.f dlassq.f dlasy2.f dorm2r.f dsteqr.f \
        dtrevc.f dtrexc.f dtrsen.f dtrsyl.f \
	ilaenv.f lsame.f lsamen.f xerbla.f xlaenv.f icmax1.f izmax1.f  \
	cgeqr2.f clacon.f clacpy.f cladiv.f clahqr.f clange.f clanhs.f \
        clarf.f  clarfg.f clarnv.f clartg.f clascl.f claset.f classq.f \
        clatrs.f cmach.f  crot.f   ctrevc.f ctrexc.f ctrsen.f ctrsyl.f \
	zgeqr2.f zlacon.f zlacpy.f zladiv.f zlahqr.f zlange.f zlanhs.f \
        zlarf.f  zlarfg.f zlarnv.f zlartg.f zlascl.f zlaset.f zlassq.f \
        zlatrs.f zmach.f  zrot.f   ztrevc.f ztrexc.f ztrsen.f ztrsyl.f \
        zunm2r.f dzsum1.f cunm2r.f scsum1.f

liblapack_arpack_la_FFLAGS = $(FFLAGS_SAV)
liblapack_arpack_la_LIBADD = liblapack_noopt.la
