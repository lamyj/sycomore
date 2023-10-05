Installing Sycomore
===================

Packaged
--------

Packaged versions of Sycomore are available on `Anaconda`_ for Linux, macOS and Windows: to install, type ``conda install -c conda-forge sycomore``.

From source
-----------

The recommended way to install from the source code is to work in a Conda environment. Git and C++14 compilers must be available; if this is not the case, they can be installed through Conda: ``conda install -c conda-forge compilers git``

Start by setting up the environment:

1. Clone the `source repository`_: ``git clone https://github.com/lamyj/sycomore.git``
2. Install the dependencies: ``conda -c conda-forge install boost boost-cpp cmake make numpy pybind11 scipy xsimd xtensor xtensor-python``

In the source directory, create a ``build`` directory and configure the build from there: ``cmake -DCMAKE_BUILD_TYPE=Release ..`` In addition to the common CMake options (e.g. *CMAKE_BUILD_TYPE* or *CMAKE_INSTALL_PREFIX*), Sycomore may be built with the following options:

- *BUILD_TESTING* controls the build of the C++ unit test executables; defaults to *ON*, i.e. the unit tests are compiled
- *BUILD_PYTHON_WRAPPERS* controls the build of the Python wrappers; defaults to *ON*, i.e. the Python wrappers are built
- *BUILD_EXAMPLES* controls the build of the C++ examples; defaults to *ON*, i.e. the C++ examples are built

Once the build is configured, run it: ``cmake --build . --target install --config Release --parallel`` 

If the install directory is not in the standard path you will need to adjust your environment to point to the install directory (*LD_LIBRARY_PATH* on Linux, *DYLD_LIBRARY_PATH* on macOS and *PATH* on Windows, *PYTHONPATH* on all platforms).

Once the compilation succeeds, the units tests can be run from the build directory:

.. code-block:: shell
  
  export SYCOMORE_TEST_DATA=$(pwd)/../tests/data
  ctest -T Test
  python3 -m unittest discover -s ../tests/python/

.. _Anaconda: https://anaconda.org/conda-forge/sycomore
.. _source repository: https://github.com/lamyj/sycomore
