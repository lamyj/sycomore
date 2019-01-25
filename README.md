# Sycomore: a configuration model library

Sycomore simulates MRI sequences using configuration models. Inspired by Carl Ganter's [CoMoTk](https://github.com/cganter/CoMoTk), Sycomore is fast and is usable both through its native C++ API and through a similar Python API. Except for the computation of derivatives, Sycomore is on-par with the features offered by CoMoTk.

## Compilation

Sycomore requires a C++11 compiler. To take full advantage of your CPU, [OpenMP](https://www.openmp.org/) is strongly recommended. If you want to validate your build of Sycomore, you should run the unit tests, which require [Boost.Test](https://www.boost.org/doc/libs/release/libs/test/). Sycomore uses [CMake](https://cmake.org/), so the simplest way to build it would be to create a `build` directory inside the source directory, run `cmake`, then run `make`:

```shell
mkdir build
cd build
cmake ..
make
```

The compilation can take advantage of a multi-core CPU either by using [make](https://www.gnu.org/software/make/) with the `-jN` flag (where `N` is the number of concurrent tasks, i.e. the number of cores) or by using [Ninja](https://ninja-build.org/).

This code benefits a lot of compilers optimization. It is recommended to build a `Release` version and use the `-ffast-math` compiler flag on gcc or LLVM. The recommended call to CMake is then:

```shell
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-march=native -ffast-math" ..
```

On macOS, the Homebrew version of OpenMP requires additional flags, namely:

```shell
cmake ⁠\
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-march=native -ffast-math" \
  -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/include" ⁠\
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/include" \
  -DOpenMP_C_LIB_NAMES="omp" -DOpenMP_CXX_LIB_NAMES="omp" \
  -DOpenMP_omp_LIBRARY=/usr/local/lib/libomp.dylib \
  ../
```

Once the compilation succeeds, the units tests can be run from the build directory:

```shell
export SYCOMORE_TEST_DATA=$(pwd)/../tests/data
ctest -T Test
```
