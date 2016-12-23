INSTALL_DIR = $(HOME)
CONFIG_FILE = $(INSTALL_DIR)/vnodelp/config/MacOSXWithProfil

include $(CONFIG_FILE)

CXXFLAGS += -I$(INSTALL_DIR)/vnodelp/include\
	-I$(INSTALL_DIR)/vnodelp/FADBAD++
LDFLAGS  += -L$(INSTALL_DIR)/vnodelp/lib 


basic: 	basic.o
	$(CXX) $(LDFLAGS) -o $@ basic.o -lvnode $(LIBS) 

clean:
	@-$(RM) *.o core.* basic