F77=$(FC)
FFLAGS=$(FCFLAGS)
FLIBS=$(FCLIBS)

AM_CPPFLAGS = -I/usr/local/include

lib_LIBRARIES = libTIDES.a

library_includedir=$(includedir)
library_include_HEADERS = minc_tides.h dp_tides.h mp_tides.h

libTIDES_a_SOURCES = commonITER.c commonITER.h \
                     doubEVENTS.c doubEVENTS.h \
                     doubITER.c doubITER.h \
                     doubNUM.c doubNUM.h doubNUMdef.h \
                     doubPOL.c doubPOL.h \
                     doubTODE.c doubTODE.h \
                     minc_tides.c minc_tides.h dp_tides.h
if HAVE_MP
  libTIDES_a_SOURCES += mpfrEVENTS.c mpfrEVENTS.h \
                        mpfrITER.c mpfrITER.h \
                        mpfrNUM.c mpfrNUM.h mpfrNUMdef.h \
                        mpfrPOL.c mpfrPOL.h \
                        mpfrTODE.c mpfrTODE.h mp_tides.h
endif

if HAVE_FORTRAN
  libTIDES_a_SOURCES += minf_tides.f
endif
