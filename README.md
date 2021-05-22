# Sycomore &mdash; an MRI simulation toolkit

[![PyPI - Wheel](https://img.shields.io/pypi/v/sycomore)](https://pypi.org/project/sycomore/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sycomore.svg)](https://anaconda.org/conda-forge/sycomore)

Sycomore is an MRI simulation toolkit providing Bloch simulation, Extended Phase Graph (EPG) (both regular and discrete, including 3D), and Configuration Model. Sycomore is a Python packge in which all computationnaly-intensive operations are run by a C++ backend, providing a very fast runtime.

Sycomore is free software, released under the [MIT license][], and its source code is available on [GitHub][].

A sample web application, using Sycomore paired with [Bokeh][] is available on [Heroku][]: it presents classical MRI experiments (RARE, RF-spoiling, slice profile with a selective sinc pulse), using the different simulation models of Sycomore.

## Installation

Packaged versions of Sycomore are available on [pypi][] and [Anaconda][] for Linux, macOS and Windows.

To install from [Anaconda][], type `conda install -c conda-forge sycomore`. To install from [pypi][], type `pip3 install sycomore` (or `pip install sycomore`). If you are installing from [pypi][] and no pre-compiled version is available for your platform, pip will try to install from the source archive.

If you need to install Sycomore from source, you will need a C++11 compiler, [CMake][], [xsimd][] and [pybind11][] to successfully build Sycomore. If you want to validate your build of Sycomore, you should run the unit tests, which require [Boost.Test][].

Additional details are provided in the [documentation][].

## Usage

The following code simulates a single repetition of a simple [RARE sequence][] using [regular EPG][] and plots the transverse magnetization of each echo.

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
- [Extended Phase Graph](https://sycomore.readthedocs.io/en/latest/epg/index.html)
- [Configuration Model](https://sycomore.readthedocs.io/en/latest/como.html)

[Anaconda]: https://anaconda.org/conda-forge/dicomifier
[Bokeh]: https://bokeh.org
[Boost.Test]: https://www.boost.org/doc/libs/release/libs/test/
[CMake]: https://cmake.org/
[documentation]: https://sycomore.readthedocs.io/en/latest/installation.html
[GitHub]: https://github.com/lamyj/sycomore/
[Heroku]: https://sycomore.herokuapp.com/
[MIT license]: https://en.wikipedia.org/wiki/MIT_License
[pybind11]: http://pybind11.readthedocs.io/
[pypi]: https://pypi.org/project/sycomore/
[RARE sequence]: https://doi.org/10.1002/mrm.1910030602
[regular EPG]: https://sycomore.readthedocs.io/en/latest/epg/regular.html
[xsimd]: https://xsimd.readthedocs.io/
