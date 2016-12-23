# set CONFIG_FILE and INSTALL_DIR

CONFIG_FILE ?= LinuxWithProfil
INSTALL_DIR ?= $(HOME)

#--- DO NOT CHANGE BELOW ---

export CONFIG_FILE
export INSTALL_DIR

all: 
	$(MAKE) -C lib
	$(MAKE) -C examples
	$(MAKE) -C tests



install: 
	@-rm -rf $(INSTALL_DIR)/vnodelp
	@-mkdir $(INSTALL_DIR)/vnodelp

	@-mkdir $(INSTALL_DIR)/vnodelp/lib
	@-cp    lib/libvnode.a $(INSTALL_DIR)/vnodelp/lib
#
	@-mkdir $(INSTALL_DIR)/vnodelp/include
	@-cp    include/*.h    $(INSTALL_DIR)/vnodelp/include
	@-cp    include/fadbad_intv.inc $(INSTALL_DIR)/vnodelp/include
#
	@-mkdir $(INSTALL_DIR)/vnodelp/config
	@-cp    config/*       $(INSTALL_DIR)/vnodelp/config
#
	@-mkdir $(INSTALL_DIR)/vnodelp/doc
	@-cp    doc/vnode.pdf  $(INSTALL_DIR)/vnodelp/doc
#	
	@-mkdir $(INSTALL_DIR)/vnodelp/FADBAD++
	@-cp -rf FADBAD++/* $(INSTALL_DIR)/vnodelp/FADBAD++
#
	@-cp    README $(INSTALL_DIR)/vnodelp/
	@-cp    LICENSE $(INSTALL_DIR)/vnodelp/
	@echo   " ***************************************"
	@echo   "  VNODE-LP installed in $(INSTALL_DIR)  "
	@echo   " ***************************************"

LIB_SRC = fadbad_ad.cc integ.cc matrixinverse.cc vtiming.cc	\
	intvfuncs.cc iho.cc					\
	fadbad_advar.cc hoe.cc					\
	debug.cc						\
	qr.cc							\
	solution.cc sortcolumns.cc			

export LIB_SRC

EXAMPLES_SRC = E1.cc scalar.cc basic.cc			\
	intermediate.cc integctrl.cc			\
	odeparam.cc integi.cc orderstudy.cc

clean:    
	@$(MAKE) -C lib clean
	@$(MAKE) -C examples clean	
	@$(MAKE) -C tests clean
