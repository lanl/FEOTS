FROM gcr.io/feots-224617/gnu-openmpi as devel 

ENV INSTALL_PATH="/usr/local/feots"
ENV FCFLAGS="-O3 -ffree-line-length-none"
ENV FC=gfortran
ENV TIME_AVG=no

RUN mkdir -p /build

COPY . /build/feots

RUN cd /build/feots && \
    FCFLAGS="-O0 -g -ffree-line-length-none -Wall -fcheck=all -ffpe-trap=invalid -fbacktrace" \
    /build/feots/configure --enable-openmp --enable-mpi || cat config.log && \
    make || cat config.log && make install

FROM gcr.io/feots-224617/gnu-openmpi 
 
COPY --from=devel /usr/local/feots /usr/local/feots

ENV PATH=/usr/local/feots/bin:$PATH
ENV FEOTS_INSTALL_PATH="/usr/local/feots"
