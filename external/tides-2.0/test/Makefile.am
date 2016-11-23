F77=$(FC)
FFLAGS=$(FCFLAGS)
FLIBS=$(FCLIBS)

LIBTIDES = $(top_builddir)/libTIDES/libTIDES.a

check_PROGRAMS = minc_kepler minc_lorenz

minc_kepler_SOURCES = test_minc_kepler.c minc_kepler.c minc_tides.h
minc_lorenz_SOURCES = test_minc_lorenz.c minc_lorenz.c minc_tides.h

minc_kepler_LDADD = $(LIBTIDES)
minc_kepler_DEPENDENCIES = $(LIBTIDES)

minc_lorenz_LDADD = $(LIBTIDES)
minc_lorenz_DEPENDENCIES = $(LIBTIDES)

if HAVE_FORTRAN
  check_PROGRAMS += minf_kepler minf_lorenz
endif

minf_kepler_SOURCES = test_minf_kepler.f minf_kepler.f
minf_lorenz_SOURCES = test_minf_lorenz.f minf_lorenz.f

minf_kepler_LDADD = $(LIBTIDES)
minf_kepler_DEPENDENCIES = $(LIBTIDES)

minf_lorenz_LDADD = $(LIBTIDES)
minf_lorenz_DEPENDENCIES = $(LIBTIDES)

check_PROGRAMS += std_kepler std_kepler_partials

std_kepler_SOURCES = test_std_kepler.c std_kepler.c std_kepler.h dp_tides.h

std_kepler_LDADD = $(LIBTIDES)
std_kepler_DEPENDENCIES = $(LIBTIDES)

std_kepler_partials_SOURCES = test_std_kepler_partials.c std_kepler_partials.c std_kepler_partials.h dp_tides.h

std_kepler_partials_LDADD = $(LIBTIDES)
std_kepler_partials_DEPENDENCIES = $(LIBTIDES)

check_PROGRAMS += std_lorenz std_lorenz_partials

std_lorenz_SOURCES = test_std_lorenz.c std_lorenz.c std_lorenz.h dp_tides.h

std_lorenz_LDADD = $(LIBTIDES)
std_lorenz_DEPENDENCIES = $(LIBTIDES)

std_lorenz_partials_SOURCES = test_std_lorenz_partials.c std_lorenz_partials.c std_lorenz_partials.h dp_tides.h

std_lorenz_partials_LDADD = $(LIBTIDES)
std_lorenz_partials_DEPENDENCIES = $(LIBTIDES)

if HAVE_MP
  check_PROGRAMS += mpfr_kepler mpfr_lorenz
else
  LIBS += -lm
endif

mpfr_kepler_SOURCES = test_mpfr_kepler.c mpfr_kepler.c mpfr_kepler.h mp_tides.h

mpfr_kepler_LDADD = $(LIBTIDES)
mpfr_kepler_DEPENDENCIES = $(LIBTIDES)

mpfr_lorenz_SOURCES = test_mpfr_lorenz.c mpfr_lorenz.c mpfr_lorenz.h mp_tides.h

mpfr_lorenz_LDADD = $(LIBTIDES)
mpfr_lorenz_DEPENDENCIES = $(LIBTIDES)

TESTS = $(check_PROGRAMS)

installcheck-local-no:
	env LIBTIDES='' LIBS='-lTIDES -lmpfr -lgmp -lm' $(MAKE) $(AM_MAKEFLAGS) -e check

installcheck-local:
	$(MAKE) $(AM_MAKEFLAGS) LIBTIDES='' LIBS='-L$(libdir) -lTIDES -lmpfr -lgmp -lm' clean $(check_PROGRAMS) check-TESTS
