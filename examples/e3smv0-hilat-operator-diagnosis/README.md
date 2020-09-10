# E3SMV0 HiLat Operator Diagnosis
This example documents the workflow for diagnosing transport operators and shows how this effort was done for the E3SMV0-HiLat parent model.

## Getting Started

### Create a FEOTS Database
To get started, you must create a FEOTS database; this is essentially a directory structure that is used to store information about the parent model, including a NetCDF file that contains the parent model mesh, impulse fields (produced in this example), impulse response functions from the parent model, and the processed transport operators. To create a FEOTS database, use the included `./mkdb.sh` script.

This script takes two arguments. The first argument is the path to where you are installing your FEOTS database(s), and the second is the name of the parent model you are working with.

For example
```
$ export FEOTS_DBROOT=/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG
$ ./mkdb $FEOTS_DBROOT
```
This will create the directory `/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG` with the following hierarchy 
```
 metadata.json
 mesh/
 ops/
 irf/response/
 irf/impulse/
```

You must save a NetCDF file, output by POP that has mesh/grid information stored within underneath the `mesh/` subdirectory and rename the file `mesh.nc`

The metadata.json file is a template that is used to describe parent model.

### Create Impulse Fields
If you have not created impulse fields associated with the parent model, you can use the provided Slurm batch script [`slurm/lanl/00-impulse.slurm`](./slurm/lanl/00-impulse.slurm).
Generation of the impulse function requires a POP NetCDF file that contains the mesh information for the parent model; the location of this file must be set in the `POP_MESHFILE` environment variable.

The impulse function and Greedy Graph coloring solution is written to your FEOTS Database `irf/impulse/` subdirectory and the extracted mesh information is written to your FEOTS Database `mesh/` subdirectory. To facilitate output from this impulse generation step to your FEOTS Database, set the environment variable `FEOTS_DBROOT` to the path of your FEOTS database (if you have not done so already).

```
$ export POP_MESHFILE=/path/to/pop/netcdf/file.nc
$ export FEOTS_DBROOT=/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG
$ sbatch --get-user-env slurm/lanl/00-impulse.slurm
```
*You may need to modify this slurm batch file for your system*


### Diagnose Transport Operators
Once you have created the impulse response fields, move them to your FEOTS database under the `irf/response/` subdirectory. Alternatively, you can create a symbolic link to those files, e.g.
```
$ ln -s /path/to/irf/output/*.nc /lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG/
```

Create a list of operators in txt file called `IRFs.txt`, e.g.
```
$ ls ${FEOTS_DBROOT}/irf/response/*.nc > IRFs.txt
```

Then, use the batch script to process the IRF files to transport and diffusion operators.
```
$ export FEOTS_DBROOT=/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-HILAT-5DAVG
$ sbatch --get-user-env slurm/lanl/01-operator-diagnosis.slurm
```

**Note**

You may need to update the number of job arrays that are launched to match the number of files that are processed. For example, use `wc -l` to report the number of operators
```
$ wc -l IRFs.txt
365
```
Then modify line 11 of `slurm/lanl/01-operator-diagnosis.slurm`
```
#SBATCH --array=1-365%2
```
The number after the `%` symbol indicates how many jobs in the job array can run at any given moment. Setting this number to a higher value can increase parallelism (Note that the operators can be processed in parallel over the IRF files).


# To DO
* Need comments on stencil type and correspondence with advection scheme in POP
* Need comments on required fields for diffusion operator generation
* Need comments on naming convention requirements for IRF fields in POP diagnostics
