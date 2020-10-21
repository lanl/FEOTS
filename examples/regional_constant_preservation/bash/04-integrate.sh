#!/bin/bash

mpirun -np 1 feots integrate ${FEOTS_FLAGS} \
                             --out ${OUTDIR} \
                             --param-file ./runtime.params
