#!/bin/bash
#
# Script for making empty FEOTS database directory structure
#
# Usage:
#
#   bash ./mkdb.sh INSTALL_PATH
#

mkdir -p ${1}/mesh
mkdir -p ${1}/irf/impulse
mkdir -p ${1}/irf/response
mkdir -p ${1}/ops

cat <<EOT >> ${2}/${1}/metadata.json
{
  "database_version": "0.0.0",
  "model_name": "",
  "model_description_url": "",
  "model_validation_url": "",
  "model_git_url": ""
  "operator_forcing_period_sec": 0,
  "n_operators": 0
}
EOT
