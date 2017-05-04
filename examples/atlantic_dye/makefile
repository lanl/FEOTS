
SHELL=/bin/bash
FEOTSDIR=$(shell pwd | sed 's/FEOTS.*/FEOTS/')
INITDIR=$(shell pwd)

GreedyColoring :
	make --directory=${FEOTSDIR}/build/ GreedyColoring
	mv ${FEOTSDIR}/build/GreedyColoring ./

GenerateMeshOnlyFile :
	make --directory=${FEOTSDIR}/build/ GenerateMeshOnlyFile
	mv ${FEOTSDIR}/build/GenerateMeshOnlyFile ./

OperatorDiagnosis :
	make --directory=${FEOTSDIR}/build/ OperatorDiagnosis
	mv ${FEOTSDIR}/build/OperatorDiagnosis ./

RegionalExtraction :
	make --directory=${FEOTSDIR}/build/ RegionalExtraction
	mv ${FEOTSDIR}/build/RegionalExtraction ./

FEOTSInitialize :
	make MODDIR=${INITDIR} --directory=${FEOTSDIR}/build/ FEOTSInitialize
	mv ${FEOTSDIR}/build/FEOTSInitialize ./

FEOTSDriver :
	make --directory=${FEOTSDIR}/build/ FEOTSDriver
	mv ${FEOTSDIR}/build/FEOTSDriver ./

clean :
	make --directory=${FEOTSDIR}/build/ clean
