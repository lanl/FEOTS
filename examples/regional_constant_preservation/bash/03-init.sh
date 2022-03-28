#!/bin/bash

set -x
mpirun -np 1 ./init --out ${OUTDIR} \
                    --regional-db ${REGIONAL_DB} \
                    --param-file ./runtime.params
