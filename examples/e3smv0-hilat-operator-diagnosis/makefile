
SHELL=/bin/bash


FC=h5pfc
OPT=-cpp -DHAVE_MPI -fopenmp -O0 -g -ffree-line-length-none -Wall -fcheck=all -ffpe-trap=invalid -fbacktrace
LIB+=$(shell nc-config --flibs) ${FEOTS_LIB}
INC+=$(shell nc-config --fflags) ${FEOTS_INCLUDE}

all :
	make init

init : FEOTSInitialize.o
	${FC} ${OPT} FEOTSInitialize.o ${INC} ${LIB} -o $@

FEOTSInitialize.o :
	${FC} ${OPT} -c FEOTSInitialize.f90 ${INC} ${LIB} -o $@	     


clean :
	rm *.o init
