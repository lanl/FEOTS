FROM gcr.io/feots-224617/gnu-openmpi as devel 

ENV INSTALL_PATH="/usr/local/feots"
ENV FCFLAGS="-O3 -ffree-line-length-none"
ENV FC=gfortran
ENV DEBUG=yes
ENV PROFILE=no
ENV NCLIB="-L/usr/local/hdf5/lib -lhdf5_fortran -lhdf5 -lhdf5_hl -L/usr/local/netcdf/lib -lnetcdff -lnetcdf -lcurl"
ENV NCINC="-I/usr/local/hdf5/include -I/usr/local/netcdf/include"
ENV OMP=yes
ENV MPI=yes
ENV TIME_AVG=no
ENV VERTICAL_DIFFUSION=yes
ENV VOLUME_CORRECTION=yes

RUN mkdir -p /build

COPY . /build/feots

RUN cd /build/feots/build && \
    make && make install

FROM gcr.io/feots-224617/gnu-openmpi 
 
COPY --from=devel /usr/local/feots /usr/local/feots

ENV FEOTS_INSTALL_PATH="/usr/local/feots"
ENV FEOTS_FCFLAGS="-O3 -ffree-line-length-none"
ENV FEOTS_FC=gfortran
ENV FEOTS_DEBUG=yes
ENV FEOTS_PROFILE=no
ENV FEOTS_OMP=yes
ENV FEOTS_MPI=yes
ENV FEOTS_TIME_AVG=no
ENV FEOTS_VERTICAL_DIFFUSION=yes
ENV FEOTS_VOLUME_CORRECTION=yes
