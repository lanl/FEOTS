
SHELL=/bin/bash


FC=h5pfc
OPT=-cpp -O3 -ffree-line-length-none
LIB+=${FEOTS_LIB} $(shell nc-config --flibs) $(shell nc-config --libs)
INC+=$(shell nc-config --cflags) $(shell nc-config --fflags) ${FEOTS_INCLUDE}


.PHONY: clean

all :
	make init
	make genmask

init : FEOTSInitialize.o
	${FC} ${OPT} FEOTSInitialize.o ${INC} ${LIB} -o $@

FEOTSInitialize.o :
	${FC} ${OPT} -c FEOTSInitialize.f90 ${INC} ${LIB} -o $@	     

genmask : GenMask.o
	${FC} ${OPT} GenMask.o ${INC} ${LIB} -o $@

GenMask.o :
	${FC} ${OPT} -c GenMask.f90 ${INC} ${LIB} -o $@	     

clean :
	rm *.o init genmask
