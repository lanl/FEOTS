# FEOTS - Regional Constant Preservation

## Description
This demo is meant to verify that constant tracer fields remain constant. This serves as a test for the volume correction
feature in FEOTS to verify that tracer is not lost through the ocean surface due to velocity field divergence/convergence
in the top-most grid cells.

The regional configuration is set up to match the [Zapiola Rise](../zapiola/README.md) demonstration. In this demonstration,
we only use one tracer field that is initialized to a value of 1 everywhere, including through the boundary conditions.
The exact solution for this configuration is that the tracer value remains at a value of 1 everywhere. Any deviation is 
a combination of numerical and finite precision errors. FEOTS was designed to preserve constant tracer fields, so that the
numerical error should be identically zero. Because of this, calculated solutions should only contain the exact solution
plus finite precision errors.

This demonstration shows how finite-precision errors from constant tracer fields can enter a solution and grow over time.


## Running this demo - "Bare Metal"
This demonstration contains two custom `.F90` files for regional mask generation and initial conditions generation.
To build these programs, you need to have the following dependencies installed 
* HDF5
* NetCDF
* FEOTS

Further, you must define the following environment variables
* `FEOTS_INSTALL` - points to the install path for feots
* `FEOTS_LIB` - points the the lib/ path for feots
* `FEOTS_INCLUDE` - points to the include/ path for feots

To build the regional mask and initial condition generation programs, you can use the provided makefile to build the programs
```
make 
```

The provided `feots.slurm` script provides a job submission script that can be used to run this demo (end-to-end).
This script

## Running this demo - Containerized
```
# Pull the latest feots docker image down from GCR
singularity pull gcr.io/feots-224617/feots:latest

# Create the regional mask file
singularity run feots.sif /opt/feots/examples/regional_constant_preservation/genmask

# Extract the regional operators from global operators (must bind-mount the FEOTS-DB to /mnt/feots-db)
singularity run feots.sif feots region-extraction --paramfile=/opt/feots/examples/regional_constant_preservation/runtime.params

# Generate Initial Conditions
mpirun -np 2 singularity run feots.sif /opt/feots/examples/regional_constant_preservation/init

# Run FEOTS in forward-integration mode
mpirun -np 2 singularity run feots.sif feots integrate --paramfile=/opt/feots/examples/regional_constant_preservation/runtime.params
```

