Installing Sycomore
===================

Packaged
--------

A packaged version of Sycomore is available on `pypi`_, either in source form or in compiled (`wheel`_) form (Linux only):

.. code-block:: shell
  
  pip3 install sycomore

Note that for platforms where no wheel version is available, you will need a C++11 compiler, `CMake`_ and `pybind11`_ to successfully build Sycomore.

From source
-----------

Sycomore requires a C++11 compiler, Python (≥ 3.5) and `pybind11`_. To take full advantage of your CPU, `OpenMP`_ is strongly recommended. If you want to validate your build of Sycomore, you should run the unit tests, which require `Boost.Test`_. Sycomore uses `CMake`_, so the simplest way to build it would be to create a *build* directory inside the source directory, run *cmake*, then run *make*:

.. code-block:: shell
  
  mkdir build
  cd build
  cmake ..
  make

To build without OpenMP support, pass *-D USE_OPENMP=OFF* to the CMake call.

The compilation can take advantage of a multi-core CPU either by using `make`_ with the *-jN* flag (where *N* is the number of concurrent tasks, i.e. the number of cores) or by using `Ninja`_. Once the compilation succeeds, the units tests can be run from the build directory:

.. code-block:: shell
  
  export SYCOMORE_TEST_DATA=$(pwd)/../tests/data
  ctest -T Test

This code benefits a lot of compilers optimization. It is recommended to build a *Release* version and use extra compiler flags to speed-up the simulations. Moreover, using OpenMP on macOS requires additional flags. Details of these flags on Linux and macOS are provided in the next sections.

Debian
......

The following call to *apt-get* will install all dependencies:

.. code-block:: shell
  
  apt-get install -y cmake g++ libboost-dev libboost-test-dev python3-pybind11

When compiling with GCC on Debian 9, *-ffast-math* enables *-ffinite-math-only*, which breaks some parts of Sycomore (the application of a time interval with gradients or the computation of isochromat with off-resonance effects). All `other optimizations turned on by -ffast-math`_ can however be used and the recommended call to CMake is:

.. code-block:: shell
  
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-march=native -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fexcess-precision=fast -D__FAST_MATH__" \
    ..

For more details on floating point math with GCC, refer to the `official documentation`_.

macOS with Homebrew
...................

The following call to `brew` will install all dependencies:

.. code-block:: shell
  
  brew install boost cmake libomp pybind11

The documentation of the `-ffast-math option in Clang`_ is rather terse, but the `source code`_ provides more details. Despite disabling non-finite maths, using *-ffast-math* does not break Sycomore. With the extra flags required to detect OpenMP, the recommended call to CMake is:

.. code-block:: shell
  
  OpenMP_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/include"
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-march=native -ffast-math" \
    -DOpenMP_C_FLAGS="${OpenMP_FLAGS}" -DOpenMP_CXX_FLAGS="${OpenMP_FLAGS}" \
    -DOpenMP_C_LIB_NAMES="omp" -DOpenMP_CXX_LIB_NAMES="omp" \
    -DOpenMP_omp_LIBRARY=/usr/local/lib/libomp.dylib \
    ../

.. _Boost.Test: https://www.boost.org/doc/libs/release/libs/test/
.. _CMake: https://cmake.org/
.. _-ffast-math option in Clang: https://clang.llvm.org/docs/UsersManual.html#cmdoption-ffast-math
.. _make: https://www.gnu.org/software/make/
.. _Ninja: https://ninja-build.org/
.. _official documentation: https://gcc.gnu.org/wiki/FloatingPointMath
.. _OpenMP: https://www.openmp.org/
.. _other optimizations turned on by -ffast-math: https://gcc.gnu.org/onlinedocs/gcc-8.2.0/gcc/Optimize-Options.html#index-ffast-math
.. _pybind11: http://pybind11.readthedocs.io/
.. _pypi: https://pypi.org/project/sycomore/
.. _source code: https://github.com/llvm-mirror/clang/blob/release_80/lib/Driver/ToolChains/Clang.cpp#L2278-L2288
.. _wheel: https://pythonwheels.com/
