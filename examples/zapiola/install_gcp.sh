#!/bin/bash

# Compiles and installs FEOTS executables on GCP

module load gcc/7.4.0 openmpi/3.1.2/gcc/7.4.0 hdf5/1.10.3/parallel/openmpi/3.1.2/gcc/7.4.0 netcdf/4.6.1/parallel/openmpi/3.1.2/gcc/7.4.0
source FEOTS_Settings.gcp

make all 
