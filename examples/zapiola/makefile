
SHELL=/bin/bash


FC=h5pfc
OPT=-cpp -O0 -g -ffree-line-length-none -Wall -fcheck=all -ffpe-trap=invalid -fbacktrace -DHAVE_MPI
#OPT=-cpp -O3 -ffree-line-length-none
LIB+=$(shell nc-config --flibs) ${FEOTS_LIB}
INC+=$(shell nc-config --fflags) ${FEOTS_INCLUDE}


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
