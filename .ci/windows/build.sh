#!/bin/sh

set -ev

export WORKSPACE=${WORKSPACE:?}
export BUILD_DIR="./build" #${WORKSPACE}/${BUILD_DIR:-build}"
export INSTALL_DIR="./install" #"${WORKSPACE}/${INSTALL_DIR:-install}"

mkdir -p ${BUILD_DIR}
mkdir -p ${INSTALL_DIR}

cd ${BUILD_DIR}

cmake \
  -DCMAKE_INSTALL_PREFIX="../${INSTALL_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -Dpybind11_DIR="C:/Libraries/pybind11/share/cmake/pybind11" \
  -DBOOST_ROOT="C:/Libraries/boost_1_67_0" \
  ${CMAKE_OPTIONS} \
  ../

cmake --build . --config Release --target install --parallel 3
