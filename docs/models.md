# Configuration Models

Once [species](species.md), [pulses](pulses.md), and [time intervals](time_intervals.md) have been defined, the simulation itself is performed by the `Model` class. A configuration model will simulate a single species, and time intervals must be defined at the creation of the model; pulses need not be defined at this time.

The creation of a configuration model requires a `Species` object, the corresponding initial magnetization, and at least one time interval:

```cpp
#include <sycomore/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/Species.h>
#include <sycomore/TimeInterval.h>
#include <sycomore/units.h>

void main()
{
    using namespace sycomore::units;
    
    sycomore::Species const water(4000_ms, 2000_ms, 2.317_um*um/s);
    sycomore::TimeInterval const half_echo(20_ms);
    sycomore::Magnetization const m0(0,0,1);
    
    sycomore::Model model(water, m0, {{"half_echo", half_echo}});
}
```

Note that the last parameter uses an initializer list to build a `std::map<std::string, sycomore::TimeInterval>`: each time interval corresponds to a dimension of the configuration model and are thus refered to by name when applying them (see following examples). Initializer lists may in fact be used for every argument; users are free to chose the syntax they find most readable:

```cpp
#include <sycomore/Model.h>
#include <sycomore/units.h>

void main()
{
    using namespace sycomore::units;
    
    sycomore::Model model(
      {4000_ms, 2000_ms, 2.317_um*um/s}, 
      {0,0,1}, 
      {{"half_echo", {500_ms}}});
}
```

## Spin echo with a hard pulse

The pulses and time intervals (with or without gradients) in a sequence can be simulated by the `apply_pulse` and `apply_time_interval` functions of the `Model` class. Both take a single argument: the `apply_pulse` function expects a `Pulse` object, while the `apply_time_interval` function expects the name of a time interval which was used in the constructor of the `Model` object. At any point of the simulation, the current magnetization can be computed from the configurations by the `isochromat` function. The following code sample shows the simulation of a simple spin echo simple, with a 90 ° excitation pulse and a 180 ° refocalization pulse, both of which are simulated by hard pulses.

```cpp
#include <iostream>

#include <sycomore/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/Species.h>
#include <sycomore/TimeInterval.h>
#include <sycomore/units.h>

void main()
{
    using namespace sycomore::units;
    
    auto const TR = 700_ms;
    auto const TE = 20_ms;
    
    sycomore::Species const water(4000_ms, 2000_ms, 2.317_um*um/s);
    
    sycomore::TimeInterval const half_echo(TE/2);
    sycomore::TimeInterval const rest(TR-TE);
    
    sycomore::Pulse const excitation(90_deg, 0_deg);
    sycomore::Pulse const refocalization(180_deg, 0_deg);
    
    sycomore::Model const model(
        water, {0,0,1}, {{"half_echo", half_echo}, {"rest", rest}});
    
    model.apply_pulse(excitation);
    model.apply_time_interval("half_echo");
    model.apply_pulse(refocalization);
    model.apply_time_interval("half_echo");
    auto const echo_signal = model.isochromat();
    
    model.apply_time_interval("rest");
    auto const end_tr_signal = model.isochromat();
    
    std::cout 
        << "Transversal magnetization at echo: " 
        << echo_signal.transversal() 
        << "; longitudinal magnetization at end of TR: " << end_tr_signal.z
        << "\n";
}
```

The `isochromat` function return a `Magnetization` object, i.e. a 3D vector with members `x`, `y` and `z` and with a `transversal()` member function which return the magnitude of the transversal magnetization. The `isochromat` function takes the following optional arguments:

- the set of configurations to include in the computation; by default all configurations are used, and an empty set (`{}`) specifies this explicitely
- the location at which to perform the computation, in order to compute e.g. slice profiles
- the relative frequency to account for off-resonance effects
