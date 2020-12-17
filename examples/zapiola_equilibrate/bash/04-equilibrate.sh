#!/bin/bash

mpirun -np 1 feots equilibrate ${FEOTS_FLAGS} \
                             --regional-db ${REGIONAL_DB} \
                             --out ${OUTDIR} \
                             --no-vertical-mixing \
                             --param-file ./runtime.params
