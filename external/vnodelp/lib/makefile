include ../config/$(CONFIG_FILE)

CXXFLAGS += -I../include 

LIB_OBJ = $(patsubst %.cc, %.o , $(LIB_SRC))

libvnode.a: ${LIB_OBJ}  
	$(AR) $(ARFLAGS) $@ $?

clean:
	@-$(RM) *.o *.h core.*

cleanall:
	@-$(RM) *.o *.cc *.h *.a core.*
