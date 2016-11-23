#
# template Makefile for development with libkhepera
#
#
.PHONY: clean


# Modify this!
KTEAM_HOME = ../../../../khepera4_development

# And maybe these
LIBKHEPERA_ROOT  = ${KTEAM_HOME}/libkhepera-1.0

TARGET_SYSTEM	= khepera-2.6

OE_HOME = /usr/local/khepera4-oetools

KHEPERA_TOOLS = ${OE_HOME}

KTEAM_KERNEL_VERSION = 2.6.24

KTEAM_KERNEL_HOME = ${KHEPERA_TOOLS}/tmp/work/overo-angstrom-linux-gnueabi/linux-omap3-${KTEAM_KERNEL_VERSION}-r97/git

# And don't touch these
ARCH = arm
CROSS_COMPILE = arm-angstrom-linux-gnueabi-
PATH := $(PATH):${KHEPERA_TOOLS}/tmp/sysroots/i686-linux/usr/armv7a/bin

CC = ${CROSS_COMPILE}gcc
CXX = ${CROSS_COMPILE}g++
LD = ${CROSS_COMPILE}ld
AR = ${CROSS_COMPILE}ar
AS = ${CROSS_COMPILE}as


INCPATH = ${LIBKHEPERA_ROOT}/build-${TARGET_SYSTEM}/include


LIBPATH = ${LIBKHEPERA_ROOT}/build-${TARGET_SYSTEM}/lib


# Pointer to the libkhepera build directory
LIBKHEPERA = ${LIBKHEPERA_ROOT}/build-${TARGET_SYSTEM}


SRCS	= $(wildcard *.c)
OBJS	= $(patsubst %.c,%.o,${SRCS})
INCS	= -I ${LIBKHEPERA}/include
LIBS	= -L ${LIBKHEPERA}/lib -lkhepera -lpthread

CFLAGS 	= -O2

#for debugging
CFLAGS 	= -g


.PHONY: all clean depend

 


SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src/\
							-I$(SCOTSROOT)/src/sparse 
						


bot_scots: bot_scots.o 
	@echo "Building $@"
	$(CXX) -o $@ $? $(LIBS) $(SCOTSINC) $(INCS) $(CFLAGS)

bot_scots-static: bot_scots.o
	@echo "Building $@"
	@$(CXX) -o $@ $? $(LIBS) -static $(INCS) $(CFLAGS) $(SCOTSINC)

all: 	${TARGET}

clean : 
	@echo "Cleaning"
	@rm -f *.o .depend ${TARGET} *~

depend:	
	@echo "Building dependencies"
	@rm -f .depend
	@touch .depend
	@makedepend ${SYS_INCLUDES} ${INCS} -Y -f .depend ${SRCS}

%.o:	%.c
	@echo "Compiling $@"
	$(CC) $(INCS) -c $(CFLAGS) $< -o $@

ifeq (.depend,$(wildcard .depend))
include .depend 
endif
