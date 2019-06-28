Sycomore -- an MRI simulation toolkit
=====================================

Sycomore is an MRI simulation toolkit providing :doc:`Bloch simulation<bloch>`, :doc:`Extended Phase Graphs (EPG)<epg/index>` (both :doc:`regular<epg/regular>` and :doc:`discrete<epg/discrete>`), and :doc:`Configuration Models<como>`. Sycomore is a Python packge in which all computationnaly-intensive operations are run by a C++ backend, providing a very fast runtime and further acceleration through `OpenMP`_.

Installation
------------

Sycomore requires a C++11 compiler, Python (≥ 3.5) and `pybind11`_. To take full advantage of your CPU, OpenMP is strongly recommended. If you want to validate your build of Sycomore, you should run the unit tests, which require `Boost.Test`_. Sycomore uses `CMake`_, so the simplest way to build it would be to create a *build* directory inside the source directory, run *cmake*, then run *make*:

.. code-block:: shell
  
  mkdir build
  cd build
  cmake ..
  make

Additional details are provided in the :doc:`documentation<installation>`.
 
Usage
-----

The following code simulates a single repetition of a simple RARE sequence and plots the transverse magnetization of each echo.

.. code-block:: python
  
  import matplotlib.pyplot
  import numpy
  import sycomore
  from sycomore.units import *

  species = sycomore.Species(1000*ms, 100*ms, 1*um**2/ms)
  TE = 4*ms
  train_length = 40

  model = sycomore.epg.Regular(species)
  data = numpy.zeros(
      train_length, dtype=[("time", sycomore.Quantity), ("signal", complex)])

  model.apply_pulse(90*deg)
  for echo in range(train_length):
     model.apply_time_interval(TE/2)
     model.apply_pulse(180*deg)
     model.apply_time_interval(TE/2)
     
     data[echo] = (((1+echo)*TE), model.echo)

  times = [x.convert_to(ms) for x in data["time"]]
  magnitude = numpy.abs(data["signal"])
  matplotlib.pyplot.plot(times, magnitude, ".", label="Simulated")
  matplotlib.pyplot.plot(
     times, [numpy.exp(-(x*species.R2).magnitude) for x in data["time"]],
     label="$T_2$ decay")

  matplotlib.pyplot.ylim(0,1)
  matplotlib.pyplot.xlabel("Time (ms)")
  matplotlib.pyplot.ylabel("Magnitude")
  matplotlib.pyplot.legend()

.. figure:: rare.png
  :alt: T2 decay in RARE
  
  T2 decay in RARE

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
