#!/bin/bash
source /opt/intel/oneapi/setvars.sh
export PATH=/opt/intel/oneapi/compiler/latest/linux/bin/intel64:/opt/intel/oneapi/compiler/latest/linux/bin:$PATH
export LD_LIBRARY_PATH=/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin

cd src/pflotran
make clean
make -j4 pflotran

# Make sure pflotran is built properly
export ARTIFACT_DIR=/tmp/test-pflotran-build
rm -Rf $ARTIFACT_DIR
mkdir -p $ARTIFACT_DIR

./pflotran -help intro
EXIT_STATUS=$?

if [ $EXIT_STATUS -eq 0 ]; then
  echo 'success' > $ARTIFACT_DIR/status
else
  echo 'failed' > $ARTIFACT_DIR/status
fi