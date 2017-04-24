# makefile
#
#
#  Directives :
#
#    

SHELL=/bin/bash
FC=gfortran
LIB=-L/home/schoonover/Software/netcdf/lib -lnetcdff -lnetcdf
INC=-I/home/schoonover/Software/netcdf/include
WFLAGS= 


ifeq (${DEBUG},yes)
  OPT=-O0 -g -ffree-line-length-none -ffpe-trap=invalid -fbacktrace
else
  OPT=-O2 -ffree-line-length-none
endif


#/////////////////////////////////////////////////////////////////////// #
FEOTSDIR=$(shell pwd | sed 's/\(FEOTS\).*/\1/g')
SRCDIR=${FEOTSDIR}/src/


all :
	make GreedyColoring
	make GenerateMeshOnlyFile
	make RegionalExtraction 
	make OperatorDiagnosis
	make FEOTSInitialize
	make FEOTSDriver
	make clean

# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Modules
ModelPrecision.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o
	${FC} ${OPT} -c ${SRCDIR}common/ConstantsDictionary.f90 -o $@

CommonRoutines.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}common/CommonRoutines.f90 ${LIB} ${INC} -o $@

BinaryIO.o : ModelPrecision.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/BinaryIO.f90 -o $@

# ---------------------------------------------------------------------- #
#                          ~/src/matrices/                               #
# ---------------------------------------------------------------------- #

CRSMatrix_Class.o : ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}matrices/CRSMatrix_Class.f90 -o $@

# ---------------------------------------------------------------------- #
#                     ~/src/solutionstorage/                             #
# ---------------------------------------------------------------------- #

TracerStorage_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o CRSMatrix_Class.o
	${FC} ${OPT} -c ${SRCDIR}solutionstorage/TracerStorage_Class.f90 -o $@

# ---------------------------------------------------------------------- #
#                              ~/src/POP/                                #
# ---------------------------------------------------------------------- #

POP_GridTypeMappings.o : ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}POP/POP_GridTypeMappings.f90 -o $@

POP_Stencil_Class.o : ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Stencil_Class.f90 -o $@

POP_AdjacencyGraph_Class.o : CommonRoutines.o POP_GridTypeMappings.o POP_Stencil_Class.o \
                             POP_Mesh_Class.o
	${FC} ${OPT} -c ${SRCDIR}POP/POP_AdjacencyGraph_Class.f90 -o $@

POP_Params_Class.o : ModelPrecision.o CommonRoutines.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Params_Class.f90 -o $@

POP_Native_Class.o : ModelPrecision.o CommonRoutines.o POP_Mesh_Class.o
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Native_Class.f90 ${LIB} ${INC} -o $@

POP_Mesh_Class.o : ModelPrecision.o ConstantsDictionary.o BinaryIO.o
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Mesh_Class.f90 ${LIB} ${INC} -o $@

POP_Regional_Class.o : ModelPrecision.o POP_Mesh_Class.o POP_AdjacencyGraph_Class.o
	${FC} ${OPT} -c ${SRCDIR}POP/POP_Regional_Class.f90 ${LIB} ${INC} -o $@

POP_FEOTS_Class.o : ModelPrecision.o ConstantsDictionary.o POP_Params_Class.o \
                    POP_Regional_Class.o TracerStorage_Class.o POP_Mesh_Class.o POP_Native_Class.o
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

RegionalExtraction : ${REOBJS}
	${FC} ${OPT} ${REOBJS} ${LIB} ${INC} -o $@

OperatorDiagnosis : ${ODOBJS}
	${FC} ${OPT} ${ODOBJS} ${LIB} ${INC} -o $@

FEOTSInitialize : ${FIOBJS}
	${FC} ${OPT} ${FIOBJS} ${LIB} ${INC} -o $@

FEOTSDriver : ${FEOBJS}
	${FC} ${OPT} ${FEOBJS} ${LIB} ${INC} -o $@

GreedyGraphColoring.o : POP_Params_Class.o POP_Mesh_Class.o POP_Stencil_Class.o \
                             POP_AdjacencyGraph_Class.o 
	${FC} ${OPT} -c ${SRCDIR}POP/programs/GreedyGraphColoring.f90 ${LIB} ${INC} -o $@

GenerateMeshOnlyFile.o : POP_Params_Class.o POP_Mesh_Class.o
	${FC} ${OPT} -c ${SRCDIR}POP/programs/GenerateMeshOnlyFile.f90 -o $@

OperatorDiagnosis.o : ModelPrecision.o CRSMatrix_Class.o POP_Params_Class.o POP_Mesh_Class.o POP_Stencil_Class.o \
                      POP_AdjacencyGraph_Class.o POP_Native_Class.o POP_GridTypeMappings.o
	${FC} ${OPT} -c ${SRCDIR}POP/programs/OperatorDiagnosis.f90 ${LIB} ${INC} -o $@

RegionalExtraction.o : CRSMatrix_Class.o POP_Mesh_Class.o POP_Stencil_Class.o POP_Regional_Class.o \
                       POP_Params_Class.o
	${FC} ${OPT} -c ${SRCDIR}POP/programs/RegionalExtraction.f90 -o $@
ifdef MODDIR
FEOTSInitialize.o : POP_Params_Class.o POP_Native_Class.o POP_FEOTS_Class.o
	${FC} ${OPT} -c ${MODDIR}/FEOTSInitialize.f90 -o $@
else
FEOTSInitialize.o : POP_Params_Class.o POP_Native_Class.o POP_FEOTS_Class.o
	${FC} ${OPT} -c ${SRCDIR}POP/programs/FEOTSInitialize.f90 -o $@
endif

FEOTSDriver.o : POP_FEOTS_Class.o
	${FC} ${OPT} -c ${SRCDIR}POP/programs/FEOTSDriver.f90 -o $@


# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod

