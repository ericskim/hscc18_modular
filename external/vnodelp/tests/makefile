# CONFIG_FILE is set in vnodelp/makefile and exported in this
# file. vnodelp/makefile calls this makefile.
include ../config/$(CONFIG_FILE)

CXXFLAGS += -I../include -DNDEBUG           
LDFLAGS  += -L../lib		

LIBS = -lvnode $(I_LIBS) $(LAPACK_LIB)		\
	$(BLAS_LIB) $(GPP_LIBS) 

TESTS = test0 test01					\
	test1 test2 test3 test4 test5 test6		\
	test_n1 test_n2 test_n3 test_n4 test_n5 test_n6	\

tests: $(TESTS)
	@./test_vnode

test0:  test0.o 
	$(CXX) $(LDFLAGS) -o $@ test0.o $(LIBS)
test01:  test01.o 
	$(CXX) $(LDFLAGS) -o $@ test01.o $(LIBS)

test1:  test1.o
	$(CXX) $(LDFLAGS) -o $@ test1.o $(LIBS)
test2:  test2.o
	$(CXX) $(LDFLAGS) -o $@ test2.o $(LIBS)
test3:  test3.o 
	$(CXX) $(LDFLAGS) -o $@ test3.o $(LIBS)
test4:  test4.o
	$(CXX) $(LDFLAGS) -o $@ test4.o $(LIBS)
test5:  test5.o
	$(CXX) $(LDFLAGS) -o $@ test5.o $(LIBS)
test6:  test6.o 
	$(CXX) $(LDFLAGS) -o $@ test6.o $(LIBS)

test_n1:  test_n1.o 
	$(CXX) $(LDFLAGS) -o $@ test_n1.o $(LIBS)
test_n2:  test_n2.o 
	$(CXX) $(LDFLAGS) -o $@ test_n2.o $(LIBS)
test_n3:  test_n3.o 
	$(CXX) $(LDFLAGS) -o $@ test_n3.o $(LIBS)
test_n4:  test_n4.o 
	$(CXX) $(LDFLAGS) -o $@ test_n4.o $(LIBS)
test_n5:  test_n5.o 
	$(CXX) $(LDFLAGS) -o $@ test_n5.o $(LIBS)
test_n6:  test_n6.o 
	$(CXX) $(LDFLAGS) -o $@ test_n6.o $(LIBS)

clean:
	@-$(RM) *.o *.out core.* $(TESTS) *.bak

