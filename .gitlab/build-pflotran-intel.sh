#!/bin/bash
source /opt/intel/oneapi/setvars.sh
export PATH=/opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64:/opt/intel/oneapi/compiler/2023.1.0/linux/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/2023.1.0/linux/compiler/lib/intel64_lin

cd src/pflotran
make clean

# gnu_code_coverage=1 gnu_runtime_checks=1 catch_warnings_as_errors=1
# intel versions
make -j4 pflotran

# Make sure pflotran is built properly
export ARTIFACT_DIR_INTEL=/tmp/test-pflotran-intel
rm -Rf $ARTIFACT_DIR_INTEL
mkdir -p $ARTIFACT_DIR_INTEL

pflotran -help intro
EXIT_STATUS=$?

if [ $EXIT_STATUS -eq 0 ]; then
  echo 'success' > $ARTIFACT_DIR_INTEL/status
else
  echo 'failed' > $ARTIFACT_DIR_INTEL/status
fi