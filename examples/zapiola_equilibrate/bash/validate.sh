#!/bin/bash
#
# Uses the job-pipeline.sh to run time integration convergence tests for this example

set -x
dts=(1080)
integrators=("euler")
GCS_DEST="gs://feots-db/simulation/zapiola_equilibrate"

export FEOTS_DBROOT=/home/joe/apps/feots_output/E3SMV0-HILAT-5DAVG
export OUTDIR_ROOT=/home/joe/apps/feots_output/zapiola_equilibrate
export REGIONAL_DB=/home/joe/apps/feots_output/argentine_basin/
export FEOTS_FLAGS=''

k=0
for dt in ${dts[@]}; do
  for integrator in ${integrators[@]}; do

    sed -i -e "s/timeStepScheme.*/timeStepScheme = '$integrator',/" ./runtime.params
    export OUTDIR="${OUTDIR_ROOT}/${integrator}-${dt}"
    bash ./bash/job-pipeline.sh > ./feots.logs
#    mv ./feots.logs ${OUTDIR}

    k+=1
#    gsutil cp -r ${OUTDIR} ${GCS_DEST}
#rm -rf ${OUTDIR}

  done
done
