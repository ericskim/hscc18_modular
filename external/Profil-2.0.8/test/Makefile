BASEDIR = ..
include $(BASEDIR)/Host.cfg

INSTALLDIR  = .
OBJECTS     = TBias.o TBiasF.o
EXAMPLELIBS = Profil Bias
EXAMPLES    = TBias TBiasF

ifeq ($(SPARSE), y)
OBJECTS     += TBiasS.o TProfilS.o
EXAMPLELIBS += lr m amd
EXAMPLES    += TBiasS TProfilS
endif

MOREINCLUDES = $(ARCHDIR)
MORELIBDIRS  = $(AMDLIB)
MORECLEANS   = *.log

include $(BASEDIR)/src/Rules.mak

TBias.o:    $(BASEDIR)/lib/libBias.a
TBiasF.o:   $(BASEDIR)/lib/libBias.a
TBiasS.o:   $(BASEDIR)/lib/libBias.a
TProfilS.o: $(BASEDIR)/lib/libBias.a $(BASEDIR)/lib/libProfil.a

check: $(addsuffix $(TESTEXT),$(EXAMPLES))
	@./check.sh $^
