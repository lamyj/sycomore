Sycomore -- an MRI simulation toolkit
=====================================

Sycomore is an MRI simulation toolkit providing :doc:`Bloch simulation<bloch>`, :doc:`Extended Phase Graphs (EPG)<epg/index>` (both :doc:`regular<epg/regular>` and :doc:`discrete<epg/discrete>`), and :doc:`Configuration Models<como>`. Sycomore is a Python packge in which all computationnaly-intensive operations are run by a C++ backend, providing a very fast runtime and further acceleration through `OpenMP`_.

Installation
------------

Sycomore requires a C++11 compiler, Python (≥ 3.5) and `pybind11`_. To take full advantage of your CPU, OpenMP is strongly recommended. If you want to validate your build of Sycomore, you should run the unit tests, which require `Boost.Test`_. Sycomore uses `CMake`_, so the simplest way to build it would be to create a `build` directory inside the source directory, run *cmake*, then run *make*:

.. code-block:: shell
  
  mkdir build
  cd build
  cmake ..
  make

Additional details are provided in the :doc:`documentation<installation>`.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :titlesonly:

   installation.rst
   common_features.rst
   bloch.rst
   epg/index.rst
   como.rst
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. _Boost.Test: https://www.boost.org/doc/libs/release/libs/test/
.. _CMake: https://cmake.org/
.. _OpenMP: https://www.openmp.org/
.. _pybind11: http://pybind11.readthedocs.io/
