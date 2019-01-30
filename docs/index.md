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

The following code simulates a single repetition of a simple RARE sequence and prints the time and transverse magnetization of each echo.

```cpp
#include <iostream>

#include <sycomore/magnetization.h>
#include <sycomore/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/Species.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;

    sycomore::Species const gray_matter(1000_ms, 100_ms, 1_um*um/ms);
    auto const TR = 500_ms;
    auto const TE = 4_ms;
    auto const train_length = 40;
    sycomore::TimeInterval const half_echo(TE/2);
    sycomore::TimeInterval const rest(TR-train_length*TE);
    sycomore::Model model(
        gray_matter, {0,0,1}, {{"half_echo", half_echo}, {"rest", rest}});

    std::vector<std::pair<sycomore::Real, sycomore::Real>> signal;

    sycomore::Real time=0;
    model.apply_pulse({90_deg, 0_deg});
    for(int echo=0; echo<train_length; ++echo)
    {
        model.apply_time_interval("half_echo");
        time += half_echo.duration;

        model.apply_pulse({180_deg, 0_deg});

        model.apply_time_interval("half_echo");
        time += half_echo.duration;

        auto const m = model.isochromat();
        signal.emplace_back(time, m.transversal());
    }
    model.apply_time_interval("rest");
    time += rest.duration;


    for(auto it = signal.begin(); it!=signal.end(); ++it)
    {
        std::cout << it->first << " " << it->second << "\n";
    }
}
```

We can then witness the expected T<sub>2</sub> decay, here with Python and Matplotlib:

```python
import sys

import matplotlib.pyplot
import numpy

data = numpy.loadtxt("signal.txt")
simulated, = matplotlib.pyplot.plot(
    data[:,0], data[:,1], ".", label="Simulated")
expected, = matplotlib.pyplot.plot(
    data[:,0], numpy.exp(-data[:,0]/0.1), label="T_2 decay")
matplotlib.pyplot.ylim(0,1)
matplotlib.pyplot.xlabel("Time (s)")
matplotlib.pyplot.ylabel("Signal")
matplotlib.pyplot.legend(handles=[simulated, expected])
matplotlib.pyplot.show()
```

![T2 decay in RARE](rare.png "T2 decay in RARE")

The features and data structures are described in the documentation:

- [Units](units.md)
- [Species](species.md)
- [RF pulses](pulses.md)
- [Time intervals](time_intervals.md)
- [Configuration models](models.md)
