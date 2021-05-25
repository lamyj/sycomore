#!/bin/sh

set -eo pipefail

# Install xsimd
curl -LO https://github.com/xtensor-stack/xsimd/archive/refs/tags/7.5.0.tar.gz
tar -x -z -f 7.5.0.tar.gz
cd xsimd-7.5.0
cmake -DCMAKE_INSTALL_PREFIX=/usr/local .
sudo make install
rm -rf xsimd-7.5.0 7.5.0.tar.gz

# Install pybind11
curl -LO https://github.com/pybind/pybind11/archive/refs/tags/v2.6.2.tar.gz
tar -x -z -f v2.6.2.tar.gz
cd pybind11-2.6.2
cmake -DPYBIND11_TEST=OFF -DCMAKE_INSTALL_PREFIX=/usr/local .
make install
rm -rf pybind11-2.6.2 v2.6.2.tar.gz
