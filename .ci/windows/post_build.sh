#!/bin/sh

set -ev

export WORKSPACE=${WORKSPACE:?}
export BUILD_DIR="${WORKSPACE}/${BUILD_DIR:-build}"
export INSTALL_DIR="${WORKSPACE}/${INSTALL_DIR:-install}"

PYTHON=$(awk -F= '$0 ~ /^PYTHON_EXECUTABLE:/ { print $2 }' ${BUILD_DIR}/CMakeCache.txt)

# export DYLD_LIBRARY_PATH=${INSTALL_DIR}/lib
# export PYTHONPATH=${INSTALL_DIR}/$(${PYTHON} -c "from distutils.sysconfig import *; print(get_python_lib(True, prefix=''))")

export SYCOMORE_TEST_DATA=${WORKSPACE}/tests/data

cd "${BUILD_DIR}"
ctest -T Test
${PYTHON} -m unittest discover -s ${WORKSPACE}/tests/python/
