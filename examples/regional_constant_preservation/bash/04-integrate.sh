#!/bin/bash

set -x
mpirun -np 1 feots integrate ${FEOTS_FLAGS} \
                             --out ${OUTDIR} \
                             --regional-db ${REGIONAL_DB} \
                             --param-file ./runtime.params
