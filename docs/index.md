# Sycomore &mdash; an MRI simulation toolkit

Sycomore is an MRI simulation toolkit using the configuration model formalism. Inspired by Carl Ganter's [CoMoTk](https://github.com/cganter/CoMoTk), Sycomore is fast and is usable both through its native C++ API and through a similar Python API. Sycomore includes:

- Radiofrequency pulses
- Gradients
- Off-resonance effects
- Susceptiblity effects
- Slice profiles

## Installation

Sycomore requires a C++11 compiler. To take full advantage of your CPU, [OpenMP](https://www.openmp.org/) is strongly recommended. If you want to validate your build of Sycomore, you should run the unit tests, which require [Boost.Test](https://www.boost.org/doc/libs/release/libs/test/). Sycomore uses [CMake](https://cmake.org/), so the simplest way to build it would be to create a `build` directory inside the source directory, run `cmake`, then run `make`:

```sh
mkdir build
cd build
cmake ..
make
```

Additional details are provided in the [documentation](installation.md).

## Usage

The features and data structures are described in the documentation:

- [Units](units.md)
- [Species](species.md)
- [RF pulses](pulses.md)
- [Time intervals](time_intervals.md)
- [Configuration models](models.md)
