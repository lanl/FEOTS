
SHELL=/bin/bash


FC=mpif90
OPT=-cpp -O0 -g -ffree-line-length-none -Wall -fcheck=all -ffpe-trap=invalid -fbacktrace
#OPT=-cpp -O3 -ffree-line-length-none
LIB+=$(shell nc-config --flibs) -L${HOME}/feots/lib -lfeots
INC+=$(shell nc-config --fflags) -I${HOME}/feots/include

init : FEOTSInitialize.o
	${FC} ${OPT} FEOTSInitialize.o ${INC} ${LIB} -o $@

FEOTSInitialize.o :
	${FC} ${OPT} -c FEOTSInitialize.f90 ${INC} ${LIB} -o $@	     


clean :
	rm *.o init
