
SHELL=/bin/bash
FEOTSDIR=$(shell pwd | sed 's/FEOTS.*/FEOTS/')
INITDIR=$(shell pwd)


.PHONY: GreedyColoring \
        GenerateMeshOnlyFile \
        GenMask \
        OperatorDiagnosis \
        RegionalExtraction \
        ExtractOceanState \
        FEOTSInitialize \
        FEOTSDriver

all :
	make GreedyColoring
	make GenerateMeshOnlyFile
	make GenMask
	make OperatorDiagnosis
	make RegionalExtraction
	make ExtractOceanState
	make FEOTSInitialize
	make FEOTSDriver


GreedyColoring :
	make --directory=${FEOTSDIR}/build/ GreedyColoring
	mv ${FEOTSDIR}/build/GreedyColoring ./

GenerateMeshOnlyFile :
	make --directory=${FEOTSDIR}/build/ GenerateMeshOnlyFile
	mv ${FEOTSDIR}/build/GenerateMeshOnlyFile ./

GenMask :
	make MODDIR=${INITDIR} --directory=${FEOTSDIR}/build/ GenMask
	mv ${FEOTSDIR}/build/GenMask ./

OperatorDiagnosis :
	make --directory=${FEOTSDIR}/build/ OperatorDiagnosis
	mv ${FEOTSDIR}/build/OperatorDiagnosis ./

RegionalExtraction :
	make --directory=${FEOTSDIR}/build/ RegionalExtraction
	mv ${FEOTSDIR}/build/RegionalExtraction ./

ExtractOceanState :
	make --directory=${FEOTSDIR}/build/ ExtractOceanState
	mv ${FEOTSDIR}/build/ExtractOceanState ./

FEOTSInitialize :
	make MODDIR=${INITDIR} --directory=${FEOTSDIR}/build/ FEOTSInitialize
	mv ${FEOTSDIR}/build/FEOTSInitialize ./

FEOTSDriver :
	make --directory=${FEOTSDIR}/build/ FEOTSDriver
	mv ${FEOTSDIR}/build/FEOTSDriver ./


clean :
	make --directory=${FEOTSDIR}/build/ clean
