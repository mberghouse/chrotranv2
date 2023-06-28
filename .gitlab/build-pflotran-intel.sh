#!/bin/bash
source /opt/intel/oneapi/setvars.sh
export PATH=/opt/intel/oneapi/mpi/2021.9.0/bin:/opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2023.1.0/linux/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin

cd src/pflotran
make clean
make -j4 gnu_code_coverage=1 gnu_runtime_checks=1 catch_warnings_as_errors=1 pflotran

