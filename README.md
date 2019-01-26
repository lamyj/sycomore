# Sycomore: a configuration model library

Sycomore simulates MRI sequences using configuration models. Inspired by Carl Ganter's [CoMoTk](https://github.com/cganter/CoMoTk), Sycomore is fast and is usable both through its native C++ API and through a similar Python API. Except for the computation of derivatives, Sycomore is on-par with the features offered by CoMoTk.

## Installation

Sycomore requires a C++11 compiler. To take full advantage of your CPU, [OpenMP](https://www.openmp.org/) is strongly recommended. If you want to validate your build of Sycomore, you should run the unit tests, which require [Boost.Test](https://www.boost.org/doc/libs/release/libs/test/). Sycomore uses [CMake](https://cmake.org/), so the simplest way to build it would be to create a `build` directory inside the source directory, run `cmake`, then run `make`:

```sh
mkdir build
cd build
cmake ..
make
```

Additional details are provided in the [documentation](docs/installation.md).

## Usage

The main features are described in the documentation:

- [Units](docs/units.md)
