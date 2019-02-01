# Radiofrequency Pulses


## Hard pulses

Radiofrequency pulses are described by their flip angle and phase. Instantaneous, hard pulses are described by the `Pulse` class, and a hard-pulse approximation of an arbitrary envelope can be created with the `hard_pulse_approximation` function.

A hard pulse can be created using either [units](units.md) or raw radians values for the flip angle and the phase:

```cpp
#include <sycomore/Pulse.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    // Values in base SI values
    sycomore::Pulse const foo(M_PI/4., M_PI);
    
    // Values using units
    sycomore::Pulse const bar(45_deg, 180_deg);
}
```

## Hard pulse approximations

A hard pulse approximation is created from a hard pulse model, an envelope and a support; it may optionally define a slice selection gradient, played during the time interval separating the hard pulses. The following example shows how to create a selective SINC pulse approximation, with and without slice selection gradient. In both cases, the last parameter is a name given to the time interval between the hard pulses (see the [model](models.md) documentation for more details).

```cpp
#include <sycomore/HardPulseApproximation.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    sycomore::Pulse const hard_pulse(40_deg, 0_deg);
    
    auto const pulse_duration = 1_ms;
    auto const pulse_support_size = 101;
    auto const support = sycomore::linspace(pulse_duration, pulse_support_size);
    
    int const zero_crossings = 2;
    auto const t0 = pulse_duration/(2*zero_crossings);
    auto const slice_thickness = 1_mm;
    
    sycomore::HardPulseApproximation const sinc_pulse_no_gradient(
        hard_pulse, support, sycomore::sinc_envelope(t0), "rf");
    sycomore::HardPulseApproximation const sinc_pulse_with_gradient(
        hard_pulse, support, sycomore::sinc_envelope(t0), 
        1/t0, slice_thickness, "rf");
}
```
