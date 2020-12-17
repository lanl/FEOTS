#!/bin/bash
#
# Uses the job-pipeline.sh to run time integration convergence tests for this example

set -x
dts=(1080 540 270 135)
nTimeSteps=(400 800 1600 3200)
nStepsPerDump=(400 800 1600 3200)
integrators=("euler" "ab2" "ab3")
GCS_DEST="gs://feots-db/simulation/zapiola"

export OUTDIR_ROOT=/home/joe/apps/feots_output/zapiola/
export REGIONAL_DB=/home/joe/apps/feots_output/argentine_basin/
export FEOTS_DBROOT=/home/joe/apps/feots_output/E3SMV0-HILAT-5DAVG/

export FEOTS_FLAGS=''
export OUTDIR="${OUTDIR_ROOT}/basic-test"

k=0
for dt in ${dts[@]}; do
  for integrator in ${integrators[@]}; do

    sed -i -e "s/nTracers.*/nTracers=1,/" ./runtime.params
    sed -i -e "s/timeStepScheme.*/timeStepScheme = '$integrator',/" ./runtime.params
    sed -i -e "s/dt.*/dt = $dt,/" ./runtime.params
    sed -i -e "s/nTimeSteps.*/nTimeSteps = ${nTimeSteps[$k]},/" ./runtime.params
    sed -i -e "s/nStepsPerDump.*/nStepsPerDump = ${nStepsPerDump[$k]},/" ./runtime.params
    export OUTDIR="${OUTDIR_ROOT}/${integrator}-${dt}"
    bash ./bash/job-pipeline.sh > ./feots.logs
    mv ./feots.logs ${OUTDIR}

#    gsutil cp -r ${OUTDIR} ${GCS_DEST}
#    rm -rf ${OUTDIR}

  done
  k=$((k+1))
done
