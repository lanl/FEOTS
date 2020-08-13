# E3SMV0 HiLat Operator Diagnosis
This example documents the workflow for diagnosing transport operators and shows how this effort was done for the E3SMV0-HiLat parent model.

## Getting Started

### Create a FEOTS Database
To get started, you must create a FEOTS database; this is essentially a directory structure that is used to store information about the parent model, including a NetCDF file that contains the parent model mesh, impulse fields (produced in this example), impulse response functions from the parent model, and the processed transport operators. To create a FEOTS database, use the included `./mkdb.sh` script.

This script takes two arguments. The first argument is the path to where you are installing your FEOTS database(s), and the second is the name of the parent model you are working with.

For example
```
$ ./mkdb /lustre/scratch3/turquoise/jschoonover/feots E3SMV0-5DAVG
```
This will create the directory `/lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-5DAVG` with the following hierarchy 
```
 mesh/
 ops/
 irf/response/
 irf/impulse/
```

You must save a NetCDF file, output by POP that has mesh/grid information stored within underneath the `mesh/` subdirectory and rename the file `mesh.nc`

### Create Impulse Fields
If you have not created impulse fields associated with the parent model, you can do so by running
```
sbatch slurm/lanl/00-impulse.slurm
```
*You may need to modify this slurm batch file for your system*


### Diagnose Transport Operators
Once you have created the impulse response fields, move them to your FEOTS database under the `irf/response/` subdirectory. Create a list of operators in txt file called `IRFs.txt`, e.g.
```
$ ls /lustre/scratch3/turquoise/jschoonover/feots/E3SMV0-5DAVG/irfs/response/* > IRFs.txt
```
Then, use the batch script to process the IRF files to transport and diffusion operators.
```
sbatch slurm/lanl/01-operator-diagnosis.slurm
```

For the E3SMV0-HiLat Parent model, operator diagnosis takes about 10 minutes per operator on average.

**Note**

You may need to update the number of job arrays that are launched to match the number of files that are processed. For example, use `wc -l` to report the number of operators
```
$ wc -l IRFs.txt
365
```
Then modify line 11 of `slurm/lanl/01-operator-diagnosis.slurm`
```
#SBATCH --array=1-365%1
```
The number after the `%` symbol indicates how many jobs in the job array can run at any given moment. Setting this number to a higher value can increase parallelism (Note that the operators can be processed in parallel over the IRF files).


# To DO
* Need comments on stencil type and correspondence with advection scheme in POP
* Need comments on required fields for diffusion operator generation
* Need comments on naming convention requirements for IRF fields in POP diagnostics
