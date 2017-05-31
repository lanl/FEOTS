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
FC=gfortran
LIB=-L/home/${USER}/Software/netcdf/lib -lnetcdff -lnetcdf
INC=-I/home/${USER}/Software/netcdf/include
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


#/////////////////////////////////////////////////////////////////////// #
FEOTSDIR=$(shell pwd | sed 's/\(FEOTS\).*/\1/g')
SRCDIR=${FEOTSDIR}/src/


all :
	make GreedyColoring
	make GenerateMeshOnlyFile
	make GenMask
	make RegionalExtraction 
	make OperatorDiagnosis
	make FEOTSInitialize
	make FEOTSDriver
	make clean

# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Modules
ModelPrecision.o : ${SRCDIR}common/ModelPrecision.f90
	${FC} ${OPT} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o ${SRCDIR}common/ConstantsDictionary.f90
	${FC} ${OPT} -c ${SRCDIR}common/ConstantsDictionary.f90 -o $@

CommonRoutines.o : ModelPrecision.o ConstantsDictionary.o ${SRCDIR}common/CommonRoutines.f90
	${FC} ${OPT} -c ${SRCDIR}common/CommonRoutines.f90 ${LIB} ${INC} -o $@

BinaryIO.o : ModelPrecision.o CommonRoutines.o ${SRCDIR}common/BinaryIO.f90
	${FC} ${OPT} -c ${SRCDIR}common/BinaryIO.f90 -o $@

# ---------------------------------------------------------------------- #
#                          ~/src/matrices/                               #
# ---------------------------------------------------------------------- #

CRSMatrix_Class.o : ConstantsDictionary.o CommonRoutines.o ${SRCDIR}matrices/CRSMatrix_Class.f90
	${FC} ${OPT} -c ${SRCDIR}matrices/CRSMatrix_Class.f90 -o $@

# ---------------------------------------------------------------------- #
#                     ~/src/solutionstorage/                             #
# ---------------------------------------------------------------------- #

TracerStorage_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o CRSMatrix_Class.o \
                        ${SRCDIR}solutionstorage/TracerStorage_Class.f90
	${FC} ${OPT} -c ${SRCDIR}solutionstorage/TracerStorage_Class.f90 -o $@

# ---------------------------------------------------------------------- #
#                              ~/src/POP/                                #
# ---------------------------------------------------------------------- #

POP_GridTypeMappings.o : ConstantsDictionary.o \
                         ${SRCDIR}POP/POP_GridTypeMappings.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_GridTypeMappings.f90 -o $@

POP_Stencil_Class.o : ConstantsDictionary.o \
                      ${SRCDIR}POP/POP_Stencil_Class.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Stencil_Class.f90 -o $@

POP_AdjacencyGraph_Class.o : CommonRoutines.o POP_GridTypeMappings.o POP_Stencil_Class.o \
                             POP_Mesh_Class.o ${SRCDIR}POP/POP_AdjacencyGraph_Class.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_AdjacencyGraph_Class.f90 -o $@

POP_Params_Class.o : ModelPrecision.o CommonRoutines.o ConstantsDictionary.o \
                     ${SRCDIR}POP/POP_Params_Class.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Params_Class.f90 -o $@

POP_Native_Class.o : ModelPrecision.o CommonRoutines.o POP_Mesh_Class.o \
                     ${SRCDIR}POP/POP_Native_Class.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Native_Class.f90 ${LIB} ${INC} -o $@

POP_Mesh_Class.o : ModelPrecision.o ConstantsDictionary.o BinaryIO.o \
                   ${SRCDIR}POP/POP_Mesh_Class.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Mesh_Class.f90 ${LIB} ${INC} -o $@

POP_Regional_Class.o : ModelPrecision.o POP_Mesh_Class.o POP_AdjacencyGraph_Class.o \
                      ${SRCDIR}POP/POP_Regional_Class.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Regional_Class.f90 ${LIB} ${INC} -o $@

POP_FEOTS_Class.o : ModelPrecision.o ConstantsDictionary.o POP_Params_Class.o \
                    POP_Regional_Class.o TracerStorage_Class.o POP_Mesh_Class.o POP_Native_Class.o ${SRCDIR}POP/POP_FEOTS_Class.f90
	${FC} ${OPT} -c ${SRCDIR}POP/POP_FEOTS_Class.f90 ${LIB} ${INC} -o $@

# ---------------------------------------------------------------------- #
#                        ~/src/POP/Programs/                             #
# ---------------------------------------------------------------------- #
REOBJS = ModelPrecision.o \
         ConstantsDictionary.o \
         CommonRoutines.o \
         BinaryIO.o \
         CRSMatrix_Class.o \
         POP_Mesh_Class.o \
         POP_GridTypeMappings.o \
         POP_Stencil_Class.o \
         POP_Regional_Class.o \
         POP_Params_Class.o \
         RegionalExtraction.o

EOOBJS = ModelPrecision.o \
         ConstantsDictionary.o \
         CommonRoutines.o \
         BinaryIO.o \
         CRSMatrix_Class.o \
         TracerStorage_Class.o \
         POP_Mesh_Class.o \
         POP_Regional_Class.o \
         POP_GridTypeMappings.o \
         POP_Native_Class.o \
         POP_FEOTS_Class.o \
         POP_Params_Class.o \
         ExtractOceanState.o

ODOBJS = ModelPrecision.o \
       ConstantsDictionary.o \
       CommonRoutines.o \
       CRSMatrix_Class.o \
       POP_GridTypeMappings.o \
       POP_Stencil_Class.o \
       POP_Params_Class.o \
       POP_Mesh_Class.o \
       POP_AdjacencyGraph_Class.o \
       POP_Native_Class.o \
       OperatorDiagnosis.o 

GCOBJS = ModelPrecision.o \
       ConstantsDictionary.o \
       CommonRoutines.o \
       POP_GridTypeMappings.o \
       POP_Stencil_Class.o \
       POP_Params_Class.o \
       POP_Mesh_Class.o \
       POP_AdjacencyGraph_Class.o \
       GreedyGraphColoring.o 

GMOBJS= ModelPrecision.o \
      ConstantsDictionary.o \
      CommonRoutines.o \
      POP_Params_Class.o \
      POP_Mesh_Class.o \
      GenerateMeshOnlyFile.o

GMaOBJS= ModelPrecision.o \
      ConstantsDictionary.o \
      CommonRoutines.o \
      POP_Params_Class.o \
      POP_Mesh_Class.o \
      GenMask.o

FIOBJS = ModelPrecision.o \
         ConstantsDictionary.o \
         CommonRoutines.o \
         CRSMatrix_Class.o \
         TracerStorage_Class.o \
         POP_GridTypeMappings.o \
         POP_Stencil_Class.o \
         POP_Params_Class.o \
         POP_Mesh_Class.o \
         POP_Regional_Class.o \
         POP_Native_Class.o \
         POP_FEOTS_Class.o \
         FEOTSInitialize.o 

FEOBJS = ModelPrecision.o \
         ConstantsDictionary.o \
         CommonRoutines.o \
         CRSMatrix_Class.o \
         TracerStorage_Class.o \
         POP_GridTypeMappings.o \
         POP_Stencil_Class.o \
         POP_Params_Class.o \
         POP_Mesh_Class.o \
         POP_Regional_Class.o \
         POP_Native_Class.o \
         POP_FEOTS_Class.o \
         FEOTSDriver.o 
         
GreedyColoring : ${GCOBJS}
	${FC} ${OPT} ${GCOBJS} ${LIB} ${INC} -o $@ 

GenerateMeshOnlyFile : ${GMOBJS}
	${FC} ${OPT} ${GMOBJS} ${LIB} ${INC} -o $@

GenMask : ${GMaOBJS}
	${FC} ${OPT} ${GMaOBJS} ${LIB} ${INC} -o $@

RegionalExtraction : ${REOBJS}
	${FC} ${OPT} ${REOBJS} ${LIB} ${INC} -o $@

ExtractOceanState : ${EOOBJS}
	${FC} ${OPT} ${EOOBJS} ${LIB} ${INC} -o $@

OperatorDiagnosis : ${ODOBJS}
	${FC} ${OPT} ${ODOBJS} ${LIB} ${INC} -o $@

FEOTSInitialize : ${FIOBJS}
	${FC} ${OPT} ${FIOBJS} ${LIB} ${INC} -o $@

FEOTSDriver : ${FEOBJS}
	${FC} ${OPT} ${FEOBJS} ${LIB} ${INC} -o $@

GreedyGraphColoring.o : POP_Params_Class.o POP_Mesh_Class.o POP_Stencil_Class.o \
                        POP_AdjacencyGraph_Class.o ${SRCDIR}POP/programs/GreedyGraphColoring.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/GreedyGraphColoring.f90 ${LIB} ${INC} -o $@

GenerateMeshOnlyFile.o : POP_Params_Class.o POP_Mesh_Class.o \
                        ${SRCDIR}POP/programs/GreedyGraphColoring.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/GenerateMeshOnlyFile.f90 -o $@

ifdef MODDIR
GenMask.o : POP_Params_Class.o POP_Mesh_Class.o ${MODDIR}/GenMask.f90
	${FC} ${OPT} -c ${MODDIR}/GenMask.f90 -o $@
else
GenMask.o : POP_Params_Class.o POP_Mesh_Class.o ${SRCDIR}/POP/programs/GenMask.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/GenMask.f90 -o $@
endif

OperatorDiagnosis.o : ModelPrecision.o CRSMatrix_Class.o POP_Params_Class.o POP_Mesh_Class.o POP_Stencil_Class.o \
                      POP_AdjacencyGraph_Class.o POP_Native_Class.o POP_GridTypeMappings.o ${SRCDIR}/POP/programs/OperatorDiagnosis.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/OperatorDiagnosis.f90 ${LIB} ${INC} -o $@

RegionalExtraction.o : CRSMatrix_Class.o POP_Mesh_Class.o POP_Stencil_Class.o POP_Regional_Class.o \
                       POP_Params_Class.o ${SRCDIR}/POP/programs/RegionalExtraction.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/RegionalExtraction.f90 -o $@

ExtractOceanState.o : POP_Mesh_Class.o POP_FEOTS_Class.o POP_Native_Class.o \
                      ${SRCDIR}POP/programs/ExtractOceanState.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/ExtractOceanState.f90 -o $@

ifdef MODDIR
FEOTSInitialize.o : POP_Params_Class.o POP_Native_Class.o POP_FEOTS_Class.o ${MODDIR}/FEOTSInitialize.f90
	${FC} ${OPT} -c ${MODDIR}/FEOTSInitialize.f90 -o $@
else
FEOTSInitialize.o : POP_Params_Class.o POP_Native_Class.o POP_FEOTS_Class.o ${SRCDIR}/POP/programs/FEOTSInitialize.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/FEOTSInitialize.f90 -o $@
endif

FEOTSDriver.o : POP_FEOTS_Class.o ${SRCDIR}/POP/programs/FEOTSDriver.f90
	${FC} ${OPT} -c ${SRCDIR}POP/programs/FEOTSDriver.f90 -o $@


# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod

