# Sycomore &mdash; an MRI simulation toolkit

[![Travis Build Status](https://travis-ci.com/lamyj/sycomore.svg?branch=master)](https://travis-ci.com/lamyj/sycomore)
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/3malu4i0f9ycs7ab/branch/master?svg=true)](https://ci.appveyor.com/project/lamyj/sycomore/branch/master)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/sycomore)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sycomore.svg)](https://anaconda.org/conda-forge/sycomore)

Sycomore is an MRI simulation toolkit providing Bloch simulation, Extended Phase Graphs (EPG) (both regular and discrete, including 3D), and Configuration Models. Sycomore is a Python packge in which all computationnaly-intensive operations are run by a C++ backend, providing a very fast runtime and further acceleration through [OpenMP](https://www.openmp.org/).

Sycomore is free software, released under the [MIT license](https://en.wikipedia.org/wiki/MIT_License), and its source code is available on [GitHub](https://github.com/lamyj/sycomore/).

A sample web application, using Sycomore paired with [Bokeh](https://bokeh.org) is available on [Heroku](https://sycomore.herokuapp.com/) (note that this is using the free Heroku hosting, and start-up times might be rather long; this does not impact the performance when the application is running): it presents classical MRI experiments (RARE, RF-spoiling, slice profile with a selective sinc pulse), using the different simulation models of Sycomore.

## Installation

Packaged versions of Sycomore are available on [pypi](https://pypi.org/project/sycomore/) and [Anaconda](https://www.anaconda.com/distribution/) for Linux, macOS and Windows. The following table summarizes the availability of packages according to the version of the Python interpreter.


| Operating system | conda-forge   | PyPI          |
| ---------------- | ------------- | ------------- |
| Linux            | 3.6, 3.7      | 3.5, 3.6, 3.7 |
| macOS (≥ 9)      | 3.6, 3.7      | 3.6, 3.7      |
| Windows          | not available | 3.5, 3.6, 3.7 |

To install from [Anaconda](https://www.anaconda.com/distribution/), type `conda install -c conda-forge sycomore`. To install from [pypi](https://pypi.org/project/sycomore/), type `pip3 install sycomore` (or `pip install sycomore`). If you are installing from [pypi](https://pypi.org/project/sycomore/) and no pre-compiled version is available for your platform, pip will try to install from the source archive.

If you need to install Sycomore from source, you will need a C++11 compiler, `CMake`_ and `pybind11`_ to successfully build Sycomore. If you want to validate your build of Sycomore, you should run the unit tests, which require [Boost.Test](https://www.boost.org/doc/libs/release/libs/test/).

Additional details are provided in the [documentation](https://sycomore.readthedocs.io/en/latest/installation.html).

## Usage

The following code simulates a single repetition of a simple [RARE sequence](https://doi.org/10.1002/mrm.1910030602) using [regular EPG](https://sycomore.readthedocs.io/en/latest/epg/regular.html) and plots the transverse magnetization of each echo.

```python
import matplotlib.pyplot
import numpy
import sycomore
from sycomore.units import *

species = sycomore.Species(1000*ms, 100*ms, 1*um**2/ms)
TE = 4*ms
train_length = 40

model = sycomore.epg.Regular(species)
data = numpy.zeros(train_length, dtype=[("time", sycomore.Quantity), ("signal", complex)])

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
matplotlib.pyplot.show()
```

![T2 decay in RARE](docs/rare.png "T2 decay in RARE")

The features and data structures are described in the documentation:

- [Common features](https://sycomore.readthedocs.io/en/latest/common_features.html)
- [Bloch simulation](https://sycomore.readthedocs.io/en/latest/bloch.html)
- [Extended Phase Graphs](https://sycomore.readthedocs.io/en/latest/epg/index.html)
- [Configuration Models](https://sycomore.readthedocs.io/en/latest/como.html)
