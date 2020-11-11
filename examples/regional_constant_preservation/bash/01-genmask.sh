#!/bin/bash


set -x
feots genmask --dbroot ${FEOTS_DBROOT}  \
              --out ${OUTDIR} \
              --param-file ./runtime.params

