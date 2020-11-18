#!/bin/bash

mpirun -np 6 feots equilibrate ${FEOTS_FLAGS} \
                             --regional-db ${REGIONAL_DB} \
                             --out ${OUTDIR} \
                             --param-file ./runtime.params
