# CONFIG_FILE is set in vnodelp/makefile and exported in this
# file. vnodelp/makefile calls this makefile.

include ../config/$(CONFIG_FILE)

CXXFLAGS += -I../include	# compiler flags
LDFLAGS  += -L../lib		# library flags
LIBS     = -lvnode $(I_LIBS) $(LAPACK_LIB) \
	$(BLAS_LIB) $(GPP_LIBS) # libraries

EXAMPLES = orbit vanderpol basic E1 scalar basic	\
	intermediate integctrl odeparam integi		\
	order detest_c3

examples: $(EXAMPLES)

E1:  E1.o 
	$(CXX) $(LDFLAGS) -o $@ E1.o $(LIBS)
scalar:  scalar.o 
	$(CXX) $(LDFLAGS) -o $@ scalar.o $(LIBS)
basic:  basic.o 
	$(CXX) $(LDFLAGS) -o $@ basic.o $(LIBS)
intermediate: intermediate.o 
	$(CXX) $(LDFLAGS) -o $@ intermediate.o $(LIBS)
integctrl: integctrl.o 
	$(CXX) $(LDFLAGS) -o $@ integctrl.o $(LIBS)
odeparam: odeparam.o 
	$(CXX) $(LDFLAGS) -o $@ odeparam.o $(LIBS)
integi: integi.o 
	$(CXX) $(LDFLAGS) -o $@ integi.o $(LIBS)
order: orderstudy.o 
	$(CXX) $(LDFLAGS) -o $@ orderstudy.o $(LIBS)
detest_c3: detest_c3.o 
	$(CXX) $(LDFLAGS) -o $@ detest_c3.o $(LIBS)
vanderpol: vanderpol.o 
	$(CXX) $(LDFLAGS) -o $@ vanderpol.o $(LIBS)
orbit: orbit.o 
	$(CXX) $(LDFLAGS) -o $@ orbit.o $(LIBS)
clean:
	@-$(RM) *.o *.out core.* $(EXAMPLES) 
cleanall:
	@-$(RM) *.o *.cc *.out core.* $(EXAMPLES)
