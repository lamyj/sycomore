Installing Sycomore
===================

Packaged
--------

Packaged versions of Sycomore are available on `pypi`_ and `Anaconda`_ for Linux, macOS and Windows.

To install from `Anaconda`_, type ``conda install -c conda-forge sycomore``. To install from `pypi`_, type ``pip3 install sycomore`` (or ``pip install sycomore``). If you are installing from `pypi`_ and no pre-compiled version is available for your platform, pip will try to install from the source archive; in that case you will need a C++11 compiler, `CMake`_ (≥ 3.5), `xsimd`_ and `pybind11`_ (≥ 2.2.0) to successfully build Sycomore.

As of May 2021, compatibility with Python 2 is still possible: however due to the `end of life of Python 2`_, ensuring this compatibility is not a goal of Sycomore, and no such package is distributed.

From source
-----------

Installing Sycomore from source requires a C++11 compiler, `xsimd`_, Python (≥ 3.5) and `pybind11`_. If you want to validate your build of Sycomore, you should run the unit tests, which require `Boost.Test`_. Sycomore uses `CMake`_, so the simplest way to build it would be to create a *build* directory inside the source directory, run *cmake*, then run *make*:

.. code-block:: shell
  
  mkdir build
  cd build
  cmake ..
  make

In addition to the common CMake options (e.g. *CMAKE_BUILD_TYPE* or *CMAKE_INSTALL_PREFIX*), Sycomore may be built with the following options:

- *BUILD_SHARED_LIBS* controls the generation of static libraries or share libraries; defaults to *ON*, i.e. building shared libraries
- *BUILD_TESTING* controls the build of the C++ unit test executables; defaults to *ON*, i.e. the unit tests are compiled
- *BUILD_PYTHON_WRAPPERS* controls the build of the Python wrappers; defaults to *ON*, i.e. the Python wrappers are built
- *BUILD_STANDALONE_PYTHON_WRAPPERS*, controls whether a standalone Python library is built, wihtout any pure C++ compatibility; defaults to *OFF*, i.e. both the C++ and Python libraries are built

The compilation can take advantage of a multi-core CPU either by using `make`_ with the *-jN* flag (where *N* is the number of concurrent tasks, i.e. the number of cores) or by using `Ninja`_. Once the compilation succeeds, the units tests can be run from the build directory:

.. code-block:: shell
  
  export SYCOMORE_TEST_DATA=$(pwd)/../tests/data
  ctest -T Test
  python3 -m unittest discover -s ../tests/python/

This code benefits hugely from compiler optimizations. It is recommended to build a *Release* version and use extra compiler flags to speed-up the simulations. Details of these flags on Linux and macOS are provided in the next sections.

Debian and Ubuntu
.................

The following call to *apt-get* will install most dependencies:

.. code-block:: shell
  
  apt-get install -y cmake g++ libboost-dev libboost-test-dev libpython3-dev python3-distutils python3-pybind11

As of Debian 10 (Buster) and Ubuntu 21.04 (Hirsute Hippo), xsimd is not available as a package. Its installation is however easy: after fetching and unzipping the `latest release`_, run the following commands from within the source directory:

.. code-block:: shell
  
  cmake -DCMAKE_INSTALL_PREFIX=~/local .
  cmake --build . --target install

When compiling with GCC on Debian or Ubuntu, *-ffast-math* enables *-ffinite-math-only*, which breaks some parts of the Configuration Model implementation in Sycomore (the application of a time interval with gradients or the computation of isochromat with off-resonance effects). All `other optimizations turned on by -ffast-math`_ can however be used and the recommended call to CMake is:

.. code-block:: shell
  
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-march=native -fno-math-errno -funsafe-math-optimizations -fno-rounding-math -fno-signaling-nans -fcx-limited-range -fexcess-precision=fast -D__FAST_MATH__" \
    -DPYTHON_EXECUTABLE=/usr/bin/python3 \
    ..

For more details on floating point math with GCC, refer to the `official documentation`_.

macOS with Homebrew
...................

The following call to `brew` will install all dependencies:

.. code-block:: shell
  
  brew install boost cmake pybind11 xsimd

The documentation of the `-ffast-math option in Clang`_ is rather terse, but the `source code`_ provides more details. Despite disabling non-finite maths, using *-ffast-math* does not break Sycomore. The recommended call to CMake is:

.. code-block:: shell
  
  cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-march=native -ffast-math" \
    ../

.. _Anaconda: https://www.anaconda.com/distribution/
.. _Boost.Test: https://www.boost.org/doc/libs/release/libs/test/
.. _CMake: https://cmake.org/
.. _end of life of Python 2: https://www.python.org/dev/peps/pep-0373/
.. _-ffast-math option in Clang: https://clang.llvm.org/docs/UsersManual.html#cmdoption-ffast-math
.. _latest release: https://github.com/xtensor-stack/xsimd/releases
.. _make: https://www.gnu.org/software/make/
.. _Ninja: https://ninja-build.org/
.. _official documentation: https://gcc.gnu.org/wiki/FloatingPointMath
.. _other optimizations turned on by -ffast-math: https://gcc.gnu.org/onlinedocs/gcc-8.2.0/gcc/Optimize-Options.html#index-ffast-math
.. _pybind11: http://pybind11.readthedocs.io/
.. _pypi: https://pypi.org/project/sycomore/
.. _source code: https://github.com/llvm-mirror/clang/blob/release_80/lib/Driver/ToolChains/Clang.cpp#L2278-L2288
.. _wheel: https://pythonwheels.com/
.. _xsimd: https://xsimd.readthedocs.io/
