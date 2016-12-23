#########################################################################
#
# PROFIL/BIAS Makefile Rules (GNU make)
# -------------------------------------
#
# Copyright (C) 1998 Olaf Knueppel
#               2009 Christian Keil
#
# This file is part of PROFIL/BIAS.
#
# PROFIL/BIAS is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.
#
# $Id: Rules.mak 592 2009-01-28 10:04:19Z christian $
#
#########################################################################

LIBDIR		= $(BASEDIR)/lib
INCLUDEDIR	= $(BASEDIR)/include
ARCHDIR		= $(BASEDIR)/config/$(ARCH)

SRCDIR		= $(BASEDIR)/src
INCLUDE		= . $(SRCDIR) $(SRCDIR)/Base $(SRCDIR)/Packages $(MOREINCLUDES)
LINKEXAMPLELIBS	= $(addprefix -l,$(EXAMPLELIBS))
LIBRARYBASE	= lib$(notdir $(LIBRARY)).a
LIBRARYNAME	= $(dir $(LIBRARY))$(LIBRARYBASE)

ifneq ($(strip $(LIBRARY)),)
all:			$(LIBRARYNAME)
else
all:
endif

$(LIBRARYNAME):		$(LIBRARYNAME)($(OBJECTS))
	$(RANLIB) $@

$(LIBRARYNAME)(%.o):	%.o
	$(AR) $@ $<

ifneq ($(strip $(IAS)),)
%.o:			%.c
	$(CC) $(CFLAGS) $(DEFINES) $(addprefix -I,$(INCLUDE)) $(IASFLAGS) -o $(basename $@).is -c $<
	$(IAS) -o $(basename $@).isc $(basename $@).is
	$(AS) -o $@ $(basename $@).isc
	$(RM) $(basename $@).is $(basename $@).isc

%.o:			%.C
	$(CCPLUS) $(CPLUSFLAGS) $(DEFINES) $(addprefix -I,$(INCLUDE)) $(IASFLAGS) -o $(basename $@).is -c $<
	$(IAS) -o $(basename $@).isc $(basename $@).is
	$(AS) -o $@ $(basename $@).isc
	$(RM) $(basename $@).is $(basename $@).isc
else
%.o:			%.c
	$(CC) $(CFLAGS) $(DEFINES) $(addprefix -I,$(INCLUDE)) -o $@ -c $<

%.o:			%.C
	$(CCPLUS) $(CPLUSFLAGS) $(DEFINES) $(addprefix -I,$(INCLUDE)) -o $@ -c $<
endif

%.o:			%.s
	$(AS) -o $@ $<

%$(TESTEXT):		%.o
	$(LINK) $(LINKFLAGS) $(CFLAGS) -o $@ $< -L$(LIBDIR) $(LINKEXAMPLELIBS) $(SYSLIBS)

examples:		$(addsuffix $(TESTEXT),$(EXAMPLES))

clean:
	-@$(RM) $(LIBRARYNAME) $(OBJECTS) $(addsuffix $(TESTEXT),$(EXAMPLES))

install:		all
ifneq ($(strip $(LIBRARY)),)
	-@test -d $(LIBDIR)                   || mkdir $(LIBDIR)
endif
	-@test -d $(INCLUDEDIR)               || mkdir $(INCLUDEDIR)
	-@test -d $(INCLUDEDIR)/$(INSTALLDIR) || mkdir $(INCLUDEDIR)/$(INSTALLDIR)
ifneq ($(strip $(LIBRARY)),)
	$(INSTALL) $(LIBRARYNAME) $(LIBDIR)
	$(RANLIB) $(LIBDIR)/$(LIBRARYBASE)
endif
	$(INSTALL) $(INCLUDES) $(INCLUDEDIR)/$(INSTALLDIR)

uninstall:
	-@$(RM) $(LIBDIR)/$(LIBRARYBASE)
	-@$(RM) $(addprefix $(INCLUDEDIR)/$(INSTALLDIR)/,$(INCLUDES))
	-@$(RMDIR) $(INCLUDEDIR)/$(INSTALLDIR)
	-@$(RMDIR) $(INCLUDEDIR)
	-@$(RMDIR) $(LIBDIR)

register:
	-@mkdir RCS
	@set allow_null_glob_expansion; \
	 for file in *.[cCSh] Makefile; do \
	  echo Registering $$file; \
	  ci -u -t- $$file; \
	 done

#
# Local Variables:
# mode: Makefile
# End:
#
