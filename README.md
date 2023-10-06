# Sycomore &mdash; an MRI simulation toolkit

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sycomore.svg)](https://anaconda.org/conda-forge/sycomore)

Sycomore is an MRI simulation toolkit providing [isochromat simulation](https://sycomore.readthedocs.io/en/latest/isochromat.html) and [Extended Phase Graph (EPG)](https://sycomore.readthedocs.io/en/latest/epg/index.html). Sycomore is a Python packge in which all computationnaly-intensive operations are run by a C++ backend, providing a very fast runtime.

Sycomore is free software, released under the [MIT license][], and its source code is available on [GitHub][].

## Installation

Packaged versions of Sycomore are available on [Anaconda][] for Linux, macOS and Windows.

To install from [Anaconda][], type `conda install -c conda-forge sycomore`. Additional details, including building from source, are provided in the [documentation][].

## Usage

The following code simulates a single repetition of a simple [RARE sequence][] using [regular EPG][] and plots the transverse magnetization of each echo.

```python
import numpy
import sycomore
from sycomore.units import *

species = sycomore.Species(1000*ms, 100*ms)
TE = 4*ms
train_length = 40

model = sycomore.epg.Regular(species)
signal = numpy.zeros(train_length, dtype=complex)

model.apply_pulse(90*deg)
for echo in range(train_length):
    model.apply_time_interval(TE/2)
    model.apply_pulse(180*deg)
    model.apply_time_interval(TE/2)
    
    signal[echo] = model.echo
```

![T2 decay in RARE](docs/rare.png "T2 decay in RARE")

The features and data structures are described in the documentation:

- [Installing Sycomore](https://sycomore.readthedocs.io/en/latest/installation.html)
- [Common features](https://sycomore.readthedocs.io/en/latest/common_features.html)
- [Isochromat Simulation](https://sycomore.readthedocs.io/en/latest/isochromat.html)
- [Extended Phase Graph](https://sycomore.readthedocs.io/en/latest/epg/index.html)

[Anaconda]: https://anaconda.org/conda-forge/sycomore
[documentation]: https://sycomore.readthedocs.io/en/latest/installation.html
[GitHub]: https://github.com/lamyj/sycomore/
[MIT license]: https://en.wikipedia.org/wiki/MIT_License
[RARE sequence]: https://doi.org/10.1002/mrm.1910030602
[regular EPG]: https://sycomore.readthedocs.io/en/latest/epg/regular.html
