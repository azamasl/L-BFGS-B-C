CC = arch -x86_64 gcc
#CFLAGS = -g
#CCFLAGS = -g3 -gdwarf2 -DDEBUG -g

LDFLAGS = -lm

SRC = lbfgsb.c ls/lnsr.c ls/csearch.c ls/cstep.c  print.c linpack.c miniCBLAS.c timer.c subalgorithms.c
MAIN = driver5.c

default: build run
build :
	 $(CC)  -std=c11 $(MAIN) $(SRC) -lm -o x
run : x
	./x