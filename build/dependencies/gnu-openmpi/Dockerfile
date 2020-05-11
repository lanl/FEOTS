# 
# HPC Base image
# 
# Contents:
#   HDF5 version 1.10.5
#   GNU compilers version 8.3
#   OpenMPI version 4.0.1
#   Python 2 and 3 (upstream)
# 

FROM centos:centos7 AS devel

# Python
RUN yum install -y \
        python2 \
        python3 && \
    rm -rf /var/cache/yum/*

# GNU compiler
RUN yum install -y centos-release-scl && \
    yum install -y \
        devtoolset-8-gcc \
        devtoolset-8-gcc-c++ \
        devtoolset-8-gcc-gfortran && \
    rm -rf /var/cache/yum/*
RUN update-alternatives --install /usr/bin/g++ g++ /opt/rh/devtoolset-8/root/usr/bin/g++ 30 && \
    update-alternatives --install /usr/bin/gcc gcc /opt/rh/devtoolset-8/root/usr/bin/gcc 30 && \
    update-alternatives --install /usr/bin/gcov gcov /opt/rh/devtoolset-8/root/usr/bin/gcov 30 && \
    update-alternatives --install /usr/bin/gfortran gfortran /opt/rh/devtoolset-8/root/usr/bin/gfortran 30

# OpenMPI version 4.0.1
RUN yum install -y \
        bzip2 \
        file \
        hwloc \
        make \
        numactl-devel \
        openssh-clients \
        perl \
        tar \
        wget && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://www.open-mpi.org/software/ompi/v4.0/downloads/openmpi-4.0.1.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/openmpi-4.0.1.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/openmpi-4.0.1 &&  CC=gcc CXX=g++ F77=gfortran F90=gfortran FC=gfortran ./configure --prefix=/usr/local/openmpi --disable-getpwuid --enable-orterun-prefix-by-default --without-verbs && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/openmpi-4.0.1 /var/tmp/openmpi-4.0.1.tar.bz2
ENV LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/openmpi/bin:$PATH

# HDF5 version 1.10.5
RUN yum install -y \
        bzip2 \
        file \
        make \
        wget \
        zlib-devel && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/hdf5-1.10.5.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/hdf5-1.10.5 &&  CC=mpicc CXX=mpic++ F77=mpifort F90=mpifort FC=mpifort ./configure --prefix=/usr/local/hdf5 --enable-parallel --enable-threadsafe --enable-unsupported --enable-cxx --enable-fortran && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/hdf5-1.10.5 /var/tmp/hdf5-1.10.5.tar.bz2
ENV CPATH=/usr/local/hdf5/include:$CPATH \
    HDF5_DIR=/usr/local/hdf5 \
    LD_LIBRARY_PATH=/usr/local/hdf5/lib:$LD_LIBRARY_PATH \
    LIBRARY_PATH=/usr/local/hdf5/lib:$LIBRARY_PATH \
    PATH=/usr/local/hdf5/bin:$PATH

# NetCDF version 4.7.3, NetCDF C++ version 4.3.1, NetCDF Fortran
# version 4.5.2
RUN yum install -y \
        ca-certificates \
        file \
        libcurl-devel \
        m4 \
        make \
        wget \
        zlib-devel && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://github.com/Unidata/netcdf-c/archive/v4.7.3.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/v4.7.3.tar.gz -C /var/tmp -z && \
    cd /var/tmp/netcdf-c-4.7.3 &&  CC=mpicc CXX=mpic++ F77=mpif77 F90=mpif90 FC=mpifort ./configure --prefix=/usr/local/netcdf --enable-parallel-tests && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/netcdf-c-4.7.3 /var/tmp/v4.7.3.tar.gz
ENV CPATH=/usr/local/netcdf/include:$CPATH \
    LD_LIBRARY_PATH=/usr/local/netcdf/lib:$LD_LIBRARY_PATH \
    LIBRARY_PATH=/usr/local/netcdf/lib:$LIBRARY_PATH \
    PATH=/usr/local/netcdf/bin:$PATH
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://github.com/Unidata/netcdf-cxx4/archive/v4.3.1.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/v4.3.1.tar.gz -C /var/tmp -z && \
    cd /var/tmp/netcdf-cxx4-4.3.1 &&  CC=mpicc CXX=mpic++ F77=mpif77 F90=mpif90 FC=mpifort ./configure --prefix=/usr/local/netcdf --enable-parallel-tests && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/netcdf-cxx4-4.3.1 /var/tmp/v4.3.1.tar.gz
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://github.com/Unidata/netcdf-fortran/archive/v4.5.2.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/v4.5.2.tar.gz -C /var/tmp -z && \
    cd /var/tmp/netcdf-fortran-4.5.2 && CC=mpicc CXX=mpic++ F77=mpif77 F90=mpif90 FC=mpifort ./configure --prefix=/usr/local/netcdf --enable-parallel-tests && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/netcdf-fortran-4.5.2 /var/tmp/v4.5.2.tar.gz


