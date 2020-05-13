# makefile
#
# Copyright 2017 Joseph Schoonover, Los Alamos National Laboratory
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this 
#  software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#    

SHELL=/bin/bash
WFLAGS= 


ifeq (${DEBUG},yes)
  OPT=-cpp -O0 -g -ffree-line-length-none -Wall -fcheck=all -ffpe-trap=invalid -fbacktrace
else
  OPT=-cpp -O3 -ffree-line-length-none
endif

ifeq (${OMP},yes)
  OPT+=-fopenmp -DHAVE_OPENMP
endif

ifeq (${PROFILE},yes)
  OPT+=-pg
endif

ifeq (${MPI},yes)
   FC=mpif90
   OPT+=-DHAVE_MPI
endif

ifeq (${TIME_AVG},yes)
   OPT+=-DTIME_AVG
endif

ifeq (${VERTICAL_DIFFUSION},yes)
  OPT+=-DVERTICAL_DIFFUSION
endif

ifeq (${VOLUME_CORRECTION},yes)
  OPT+=-DVOLUME_CORRECTION
endif


LIB+=${NCLIB}
INC+=${NCINC}

#/////////////////////////////////////////////////////////////////////// #
FEOTSDIR=$(shell pwd) #| sed 's/\(FEOTS\).*/\1/g')
SRCDIR=../src/


all :
	make feots

install :
	mkdir -p ${INSTALL_PATH}/bin
	mv feots ${INSTALL_PATH}/bin
	cp -r ../examples ${INSTALL_PATH}


# Modules
ModelPrecision.o : ${SRCDIR}ModelPrecision.F90
	${FC} ${OPT} -c ${SRCDIR}ModelPrecision.F90 -o $@

ConstantsDictionary.o : ModelPrecision.o ${SRCDIR}ConstantsDictionary.F90
	${FC} ${OPT} -c ${SRCDIR}ConstantsDictionary.F90 -o $@

CommonRoutines.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}CommonRoutines.F90
	${FC} ${OPT} -c ${SRCDIR}CommonRoutines.F90 ${LIB} ${INC} -o $@

BinaryIO.o : ModelPrecision.o CommonRoutines.o ${SRCDIR}BinaryIO.F90
	${FC} ${OPT} -c ${SRCDIR}BinaryIO.F90 -o $@


CRSMatrix_Class.o : ConstantsDictionary.o CommonRoutines.o ${SRCDIR}CRSMatrix_Class.F90
	${FC} ${OPT} -c ${SRCDIR}CRSMatrix_Class.F90 -o $@


TracerStorage_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o CRSMatrix_Class.o \
                        ${SRCDIR}TracerStorage_Class.F90
	${FC} ${OPT} -c ${SRCDIR}TracerStorage_Class.F90 -o $@

POP_GridTypeMappings.o : ConstantsDictionary.o \
                         ${SRCDIR}POP_GridTypeMappings.F90
	${FC} ${OPT} -c ${SRCDIR}POP_GridTypeMappings.F90 -o $@

POP_Stencil_Class.o : ConstantsDictionary.o \
                      ${SRCDIR}POP_Stencil_Class.F90
	${FC} ${OPT} -c ${SRCDIR}POP_Stencil_Class.F90 -o $@

POP_AdjacencyGraph_Class.o : CommonRoutines.o POP_GridTypeMappings.o POP_Stencil_Class.o \
                             POP_Mesh_Class.o ${SRCDIR}POP_AdjacencyGraph_Class.F90
	${FC} ${OPT} -c ${SRCDIR}POP_AdjacencyGraph_Class.F90 -o $@

POP_Params_Class.o : ModelPrecision.o CommonRoutines.o ConstantsDictionary.o \
                     ${SRCDIR}POP_Params_Class.F90
	${FC} ${OPT} -c ${SRCDIR}POP_Params_Class.F90 -o $@

POP_Native_Class.o : ModelPrecision.o CommonRoutines.o POP_Mesh_Class.o \
                     ${SRCDIR}POP_Native_Class.F90
	${FC} ${OPT} -c ${SRCDIR}POP_Native_Class.F90 ${LIB} ${INC} -o $@

POP_Mesh_Class.o : ModelPrecision.o ConstantsDictionary.o BinaryIO.o \
                   ${SRCDIR}POP_Mesh_Class.F90
	${FC} ${OPT} -c ${SRCDIR}POP_Mesh_Class.F90 ${LIB} ${INC} -o $@

POP_Regional_Class.o : ModelPrecision.o POP_Mesh_Class.o POP_AdjacencyGraph_Class.o \
                      ${SRCDIR}POP_Regional_Class.F90
	${FC} ${OPT} -c ${SRCDIR}POP_Regional_Class.F90 ${LIB} ${INC} -o $@

POP_FEOTS_Class.o : ModelPrecision.o ConstantsDictionary.o POP_Params_Class.o \
                    POP_Regional_Class.o TracerStorage_Class.o POP_Mesh_Class.o POP_Native_Class.o ${SRCDIR}POP_FEOTS_Class.F90
	${FC} ${OPT} -c ${SRCDIR}POP_FEOTS_Class.F90 ${LIB} ${INC} -o $@

FEOTS_Driver_Routines.o : ModelPrecision.o ConstantsDictionary.o POP_Params_Class.o \
                    POP_Regional_Class.o TracerStorage_Class.o POP_Mesh_Class.o POP_Native_Class.o POP_FEOTS_Class.o ${SRCDIR}FEOTS_Driver_Routines.F90
	${FC} ${OPT} -c ${SRCDIR}FEOTS_Driver_Routines.F90 ${LIB} ${INC} -o $@

FEOTS.o : FEOTS_Driver_Routines.o ${SRCDIR}FEOTS.F90
	${FC} ${OPT} -c ${SRCDIR}FEOTS.F90 ${LIB} ${INC} -o $@

feots : FEOTS.o
	${FC} ${OPT} *.o ${LIB} ${INC} -o $@
# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod

