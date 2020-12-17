#!/bin/bash

set -x
mpirun -np 1 feots genmaps --out ${OUTDIR} \
                           --regional-db ${REGIONAL_DB} \
                           --dbroot ${FEOTS_DBROOT} \
                           --param-file ./runtime.params
