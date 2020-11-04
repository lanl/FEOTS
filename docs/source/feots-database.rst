===============
Database
===============

The FEOTS ecosystem understands two types of databases

1. Global
2. Regional

Both databases are essentially a filesystem hierarchy with expected HDF5, NetCDF, and text files stored appropriately. 


Global Database
===============

A "Global" database is the database that is closely linked to the parent model. It contains the parent model mesh, impulse response functions produced by the parent model, as well as the impulse fields and the diagnose "global" transport operators. A global FEOTS database has the following directory structure::

 metadata.json
 irf/impulse/graph.h5
 irf/impulse/ImpulseFields.nc
 irf/response/IRF.*.nc
 ops/transport.*.h5
 ops/diffusion.*.h5
 mesh/mesh.nc

metadata.json
*************

The metadata.json file stores useful information that describes the parent model and some of the characteristics of the IRF fields, such as the time average period for each operator and how many operators are provided in the database. The metadata.json file has the following schema::

 {
  "model_name": string,
  "model_description_url": string,
  "model_validation_url": string
  "model_git_url": string
  "operator_tavg_period_sec": integer,
  "n_operators": integer,
  "stencil_type": string
 }



Regional Database
==================
