#!/bin/bash

# configure intel oneapi paths
source /opt/intel/oneapi/setvars.sh
export PATH=/opt/intel/oneapi/mpi/2021.9.0/bin:/opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2023.1.0/linux/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin
ifort --version
icc --version
icpc --version

# build mpich
tar -xzvf mpich-4.1.tar.gz
MPICH_DIR=/scratch/mpich-4.1
MPICH_INSTALL_DIR=$MPICH_DIR/install
cd $MPICH_DIR
# from petsc 3.19 --download-mpich=yes configure for mpich 4.1
./configure --prefix=$MPICH_INSTALL_DIR MAKE=/usr/bin/gmake --libdir=$MPICH_INSTALL_DIR/lib CC=icc CFLAGS="-g -O0 -fPIC -diag-disable=10441" AR=/usr/bin/ar ARFLAGS=cr CXX=icpc CXXFLAGS="-g -O0 -diag-disable=10441" FFLAGS="-g -O0" FC=ifort F77=ifort FCFLAGS="-g -O0" --disable-shared --with-pm=hydra --with-hwloc=embedded --enable-fast=no --enable-error-messages=all --with-device=ch3:sock --enable-g=meminit
make all; make install

# clone and build petsc
git clone https://gitlab.com/petsc/petsc.git $PETSC_DIR
cd $PETSC_DIR
git checkout $PETSC_VERSION
$MPICH_INSTALL_DIR/bin/mpicc --version
$MPICH_INSTALL_DIR/bin/mpicxx --version
$MPICH_INSTALL_DIR/bin/mpif90 --version
./configure PETSC_ARCH=petsc-arch-intel \
--with-cc=$MPICH_INSTALL_DIR/bin/mpicc \
--with-cxx=$MPICH_INSTALL_DIR/bin/mpicxx \
--with-fc=$MPICH_INSTALL_DIR/bin/mpif90 \
--COPTFLAGS='-g -O0 -diag-disable=10441' --CXXOPTFLAGS='-g -O0 -diag-disable=10441' --FOPTFLAGS='-g -O0' --with-clanguage=c --with-shared-libraries=0 --with-debugging=1 --download-hdf5=yes --with-valgrind=1 --download-parmetis=yes --download-metis=yes --download-hypre=yes --download-mpich=yes --with-c2html=0 --with-blas-lapack-dir=/opt/intel/oneapi/mkl/2023.1.0 FFLAGS='-diag-disable 5462'
make
make check
rm -Rf petsc-arch-intel/externalpackages
