###############
Getting Started
###############

**************************
Installation with Autoconf
**************************
This section of the documentation will walk you through installing FEOTS on your own system


Dependencies
============
FEOTS requires a 2008 compliant Fortran compiler in addition to the following package dependencies

*  MPI (OpenMPI, MPICH, or MVAPICH)
*  hdf5 parallel (v1.10.x or greater)
*  NetCDF-C
*  NetCDF-Fortran

Build
==============
FEOTS comes with an autoconf build system that will autodetect dependencies by examining your PATH and LD_LIBRARY_PATH environment variables. To build and install FEOTS::

 autoreconf --install
 ./configure --prefix=/path/to/install
 make
 make install


******
Docker
******


