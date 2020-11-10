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

cat <<EOT >> ${1}/metadata.json
{
  "database_version": "0.0.1",
  "attribution": {
    "creator_names": [],
    "creator_emails": [],
    "creator_organizations": []
  },
  "hosting_url":[],
  "parent_model":{
    "name": "",
    "description_url": "",
    "validation_url": "",
    "git_url": "",
    "nominal_resolution_km": 0.0,
    "advection_scheme": "LaxWendroff27"
  },
  "operators": {
    "period_sec": 0.0,
    "n_operators": 0
  }
}
EOT
