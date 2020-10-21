#!/bin/bash

feots region-extraction --dbroot ${FEOTS_DBROOT} \
                        --out ${OUTDIR} \
                        --oplevel ${SLURM_ARRAY_TASK_ID} \
                        --param-file ./runtime.params
