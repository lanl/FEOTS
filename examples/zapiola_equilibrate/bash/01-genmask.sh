#!/bin/bash

mkdir -p ${OUTDIR}
IRF=$(sed -n 1p ${IRF_FILE})
./genmask --dbroot ${FEOTS_DBROOT}  \
          --irf ${IRF} \
          --out ${OUTDIR} \
          --param-file ./runtime.params

