#!/bin/bash

mpirun -np 1 feots integrate ${FEOTS_FLAGS} \
                             --regional-db ${REGIONAL_DB} \
                             --out ${OUTDIR} \
                             --param-file ./runtime.params
