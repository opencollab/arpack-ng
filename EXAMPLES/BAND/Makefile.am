LDADD = $(top_builddir)/SRC/libarpack$(LIBSUFFIX).la $(LAPACK_LIBS) $(BLAS_LIBS)

SSBDR = ssbdr1 ssbdr2 ssbdr3 ssbdr4 ssbdr5 ssbdr6
DSBDR = dsbdr1 dsbdr2 dsbdr3 dsbdr4 dsbdr5 dsbdr6
SNBDR = snbdr1 snbdr2 snbdr3 snbdr4 snbdr5 snbdr6
DNBDR = dnbdr1 dnbdr2 dnbdr3 dnbdr4 dnbdr5 dnbdr6
CNBDR = cnbdr1 cnbdr2 cnbdr3 cnbdr4
ZNBDR = znbdr1 znbdr2 znbdr3 znbdr4
BAND = $(SSBDR) $(DSBDR) $(SNBDR) $(DNBDR) $(CNBDR) $(ZNBDR)

check_PROGRAMS = $(BAND)
TESTS = $(check_PROGRAMS)

EXTRA_DIST = README

# Simple symmetric problem using BAND solver (single precision)
ssbdr1_SOURCES = ssbdr1.f ssband.f
ssbdr2_SOURCES = ssbdr2.f ssband.f
ssbdr3_SOURCES = ssbdr3.f ssband.f
ssbdr4_SOURCES = ssbdr4.f ssband.f
ssbdr5_SOURCES = ssbdr5.f ssband.f
ssbdr6_SOURCES = ssbdr6.f ssband.f

# Simple symmetric problem using BAND solver (double precision)
dsbdr1_SOURCES = dsbdr1.f dsband.f
dsbdr2_SOURCES = dsbdr2.f dsband.f
dsbdr3_SOURCES = dsbdr3.f dsband.f
dsbdr4_SOURCES = dsbdr4.f dsband.f
dsbdr5_SOURCES = dsbdr5.f dsband.f
dsbdr6_SOURCES = dsbdr6.f dsband.f

# Simple nonsymmetric problem using BAND solver (single precision)
snbdr1_SOURCES = snbdr1.f snband.f
snbdr2_SOURCES = snbdr2.f snband.f
snbdr3_SOURCES = snbdr3.f snband.f
snbdr4_SOURCES = snbdr4.f snband.f
snbdr5_SOURCES = snbdr5.f snband.f
snbdr6_SOURCES = snbdr6.f snband.f

# Simple nonsymmetric problem using BAND solver (double precision)
dnbdr1_SOURCES = dnbdr1.f dnband.f
dnbdr2_SOURCES = dnbdr2.f dnband.f
dnbdr3_SOURCES = dnbdr3.f dnband.f
dnbdr4_SOURCES = dnbdr4.f dnband.f
dnbdr5_SOURCES = dnbdr5.f dnband.f
dnbdr6_SOURCES = dnbdr6.f dnband.f

# Complex problem using BAND solver (single precision complex)
cnbdr1_SOURCES = cnbdr1.f cnband.f
cnbdr2_SOURCES = cnbdr2.f cnband.f
cnbdr3_SOURCES = cnbdr3.f cnband.f
cnbdr4_SOURCES = cnbdr4.f cnband.f

# Complex problem using BAND solver (double precision complex)
znbdr1_SOURCES = znbdr1.f znband.f
znbdr2_SOURCES = znbdr2.f znband.f
znbdr3_SOURCES = znbdr3.f znband.f
znbdr4_SOURCES = znbdr4.f znband.f
