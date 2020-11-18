#!/bin/bash

mpirun -np 1 feots genmaps --out ${OUTDIR} \
                           --regional-db ${REGIONAL_DB} \
                           --dbroot ${FEOTS_DBROOT} \
                           --param-file ./runtime.params

mpirun -np 6 ./init --out ${OUTDIR} \
                    --regional-db ${REGIONAL_DB} \
                    --dbroot ${FEOTS_DBROOT} \
                    --param-file ./runtime.params
