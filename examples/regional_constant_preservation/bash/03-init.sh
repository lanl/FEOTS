#!/bin/bash

mpirun -np 1 ./init --out ${OUTDIR} \
                    --param-file ./runtime.params
