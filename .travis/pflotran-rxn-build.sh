#!/bin/sh

export PETSC_DIR=$PWD/petsc; 
export PETSC_ARCH=petsc-arch

if [ $CMAKE_BUILD -eq 0 ]; then
  cd src/pflotran;
  make pflotran_rxn;
fi

