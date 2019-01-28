# Radiofrequency Pulses

Radiofrequency pulses are described by their flip angle and phase. Instantaneous, hard pulses are described by the `Pulse` class, and a hard-pulse approximation of an arbitrary envelope can be created with the `hard_pulse_approximation` function.

A hard pulse can be created using either [units](units.md) or raw radians values for the flip angle and the phase:

```cpp
#include <sycomore/Pulse.h>
#include <sycomore/units.h>

void main()
{
    using namespace sycomore::units;
    
    // Values in base SI values
    sycomore::Pulse const foo(M_PI/4., M_PI);
    
    // Values using units
    sycomore::Pulse const bar(45_deg, 180_deg)
}
```

A hard pulse approximation is created from a hard pulse model, an envelope and a support. The following example shows how to create a selective SINC pulse approximation. Note that at this points, neither timing nor gradient values are specified.

```cpp
#include <sycomore/Pulse.h>
#include <sycomore/units.h>

void main()
{
    using namespace sycomore::units;
    
    sycomore::Pulse const hard_pulse(40_deg, 0_deg);
    auto const sinc_approximation = sycomore::hard_pulse_approximation(
        hard_pulse, 
        [](sycomore::Real x) { return x==0?1:std::sin(x*M_PI)/(x*M_PI); },
        sycomore::linspace(-2., 2., 101);
}
```
