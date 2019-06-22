# Configuration Models

Once [species](species.md), [pulses](pulses.md), and [time intervals](time_intervals.md) have been defined, the simulation itself is performed by the `Model` class. A configuration model will simulate a single species, and time intervals must be defined at the creation of the model; pulses need not be defined at this time.

The creation of a configuration model requires a `Species` object, the corresponding initial magnetization, and at least one time interval:

```cpp
#include <sycomore/como/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/Species.h>
#include <sycomore/TimeInterval.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    sycomore::Species const water(4000_ms, 2000_ms, 2.317_um*um/s);
    sycomore::TimeInterval const half_echo(20_ms);
    sycomore::Magnetization const m0{0,0,1};
    
    sycomore::como::Model model(water, m0, {{"half_echo", half_echo}});
}
```

Note that the last parameter uses an initializer list to build a `std::map<std::string, sycomore::TimeInterval>`: each time interval corresponds to a dimension of the configuration model and are thus refered to by name when applying them (see following examples). Initializer lists may in fact be used for every argument; users are free to chose the syntax they find most readable:

```cpp
#include <sycomore/como/Model.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    sycomore::como::Model model(
      {4000_ms, 2000_ms, 2.317_um*um/s}, 
      {0,0,1}, 
      {{"half_echo", {500_ms}}});
}
```

## Spin echo with a hard pulse

The pulses and time intervals (with or without gradients) in a sequence can be simulated by the `apply_pulse` and `apply_time_interval` functions of the `Model` class. Both take a single argument: the `apply_pulse` function expects a `Pulse` object, while the `apply_time_interval` function expects the name of a time interval which was used in the constructor of the `Model` object. At any point of the simulation, the current magnetization can be computed from the configurations by the `isochromat` function. The following code sample shows the simulation of a simple spin echo simple, with a 90 ° excitation pulse and a 180 ° refocalization pulse, both of which are simulated by hard pulses.

```cpp
#include <iostream>

#include <sycomore/como/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/Species.h>
#include <sycomore/TimeInterval.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    auto const TR = 700_ms;
    auto const TE = 20_ms;
    
    sycomore::Species const water(4000_ms, 2000_ms, 2.317_um*um/s);
    
    sycomore::TimeInterval const half_echo(TE/2);
    sycomore::TimeInterval const idle(TR-TE);
    
    sycomore::Pulse const excitation(90_deg, 0_deg);
    sycomore::Pulse const refocalization(180_deg, 0_deg);
    
    sycomore::como::Model model(
        water, {0,0,1}, {{"half_echo", half_echo}, {"idle", idle}});
    
    model.apply_pulse(excitation);
    model.apply_time_interval("half_echo");
    model.apply_pulse(refocalization);
    model.apply_time_interval("half_echo");
    auto const echo_signal = model.isochromat();
    
    model.apply_time_interval("idle");
    auto const end_tr_signal = model.isochromat();
    
    std::cout 
        << "Transversal magnetization at echo: " 
        << sycomore::transversal(echo_signal)
        << "; longitudinal magnetization at end of TR: " << end_tr_signal[2]
        << "\n";
}
```

The `isochromat` function return a `Magnetization` object, i.e. a 3D vector with members `x`, `y` and `z` and with a `transversal()` member function which return the magnitude of the transversal magnetization. The `isochromat` function takes the following optional arguments:

- the set of configurations to include in the computation; by default all configurations are used, and an empty set (`{}`) specifies this explicitely
- the location at which to perform the computation, in order to compute e.g. slice profiles
- the relative frequency to account for off-resonance effects

## Spin echo with a slice-selective sinc pulse

The previous example can be easily adapted to use a slice-selective sinc pulse with a slice-rephasing gradient.

```cpp
#include <sycomore/como/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/Species.h>
#include <sycomore/TimeInterval.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    auto const TR = 700_ms;
    auto const TE = 20_ms;
    
    sycomore::Species const water(4000_ms, 2000_ms, 2.317_um*um/s);
    
    // Hard pulse model
    sycomore::Pulse const pulse(90_deg, 0_rad);
    
    // Approximation parameters
    auto const pulse_duration=1_ms;
    int const pulse_support_size = 101;
    int const zero_crossings = 2;
    auto const slice_thickness=1_mm;
    
    auto const t0 = pulse_duration/(2*zero_crossings);
    sycomore::HardPulseApproximation const excitation(
        pulse,
        sycomore::linspace(pulse_duration, pulse_support_size),
        sycomore::sinc_envelope(t0), 1/t0, slice_thickness, "rf");
    
    sycomore::Pulse const refocalization(180_deg, 0_rad);
    
    // Slice-rephasing gradient lobe: without ramps on the gradient, and with a
    // sinc pulse, the area of the rephasing gradient is -1/2 the area of the 
    // slice selection gradient.
    sycomore::como::Model model(
        water, {0,0,1}, {
            {excitation.get_name(), excitation.get_time_interval()},
            {"slice_rephasing", {
                (TE-pulse_duration)/2., -excitation.get_gradient_moment()/2}},
            {"half_echo", {(TE-pulse_duration)/2.}},
            {"idle", {TR-TE-pulse_duration}}
    });
    
    model.apply_pulse(excitation);
    model.apply_time_interval("slice_rephasing");
    model.apply_pulse(refocalization);
    model.apply_time_interval("half_echo");
    auto const echo_signal = model.isochromat();
    model.apply_time_interval("idle");
    auto const end_tr_signal = model.isochromat();
}
```

With the exception of definition of the hard pulse approximation, this example is very similar to the previous one: objects of type `HardPulseApproximation` are applied in the same way as objects of type `Pulse`. The last parameter of the `HardPulseApproximation` constructor is the name of the time interval which will be applied during the hard pulses when calling `apply_pulse` with a `HardPulseApproximation` object.

The previous-to-last two parameters of the `HardPulseApproximation` constructor (bandwidth and slice thickness) are optional: if they are not specified, the pulse is played without a slice selection gradient.
