Sycomore -- an MRI simulation toolkit
=====================================

Sycomore is an MRI simulation toolkit providing :doc:`isochromat simulation<isochromat>` and :doc:`Extended Phase Graph (EPG)<epg/index>` (both :doc:`regular<epg/regular>` and :doc:`discrete<epg/discrete>`, including :doc:`3D<epg/discrete_3d>`). Sycomore is a Python package in which all computationally-intensive operations are run by a C++ backend, providing a very fast runtime.

Sycomore is free software, released under the `MIT license`_, and its source code is available on `GitHub`_.

Installation
------------

Packaged versions of Sycomore are available on `Anaconda`_ for Linux, macOS and Windows.

To install from `Anaconda`_, type ``conda install -c conda-forge sycomore``. Additional details, including building from source, are provided in the :doc:`documentation<installation>`.
 
Usage
-----

The following code simulates a single repetition of a simple `RARE sequence`_ with :doc:`regular EPG<epg/regular>` and plots the transverse magnetization of each echo.

.. tab:: Python
    
    .. literalinclude:: ../examples/rare.py
      :lines: 1-18

.. tab:: C++
    
    .. literalinclude:: ../examples/rare.cpp
        :language: cpp

.. figure:: ./rare.png

   T2 decay in RARE

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :titlesonly:

   installation.rst
   common_features.rst
   isochromat.rst
   epg/index.rst
   api/index.rst

Indices and tables
==================

* :ref:`genindex`

.. _Anaconda: https://anaconda.org/conda-forge/sycomore
.. _GitHub: https://github.com/lamyj/sycomore/
.. _MIT license: https://en.wikipedia.org/wiki/MIT_License
.. _RARE sequence: https://doi.org/10.1002/mrm.1910030602
