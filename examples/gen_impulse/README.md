# Impulse Generation Demo



## Getting Started (with Docker)
* Install a POP Mesh file and runtime.params parameter file in `/path/to/input`
* Designate an output directory `/path/to/output/`

```
docker run --mount type=bind,source=/path/to/input,destination=/input \
           --volume /path/to/output:/output \
           gcr.io/feots-224617/feots:gnu-openmpi \
           feots impulse --param-file "/input/runtime.params"
```
