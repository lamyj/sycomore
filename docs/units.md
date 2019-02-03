# Units

MRI simulations deal with various quantities: times, frequencies and angular frequencies, magnetic field strength, gradient moments, etc. All those quantities use different usual prefixes, but not always the same: relaxation time and sequence timing are usually expressed in milliseconds, while the gyromagnetic ratio is expressed in _MHz/T_ (i.e. inverse of microseconds) gradients amplitudes are expressed in _mT/m_ and slice thickness in millimeters (on clinical scanners) or micrometers (on pre-clinical scanners). This wealth of units makes it very easy to get quantities wrong by a factor of 1000 or more.

Sycomore provides a unit system so that users do not have to convert their quantities to a specific unit. In C++, units may be declared either by using a custom suffix (e.g. `500_ms`) or by multiplying or dividing by the unit name (e.g. `500*ms`). Those two syntaxes can be mixed in order to use more complex units (e.g. `267.522_MHz/T`). Unit objects follow the usual arithmetic rules, and all SI [base units](https://en.wikipedia.org/wiki/SI_base_unit), [derived units](https://en.wikipedia.org/wiki/SI_derived_unit) and [prefixes](https://en.wikipedia.org/wiki/Metric_prefix). Moreover, Sycomore defines MRI-specific units, such as the diffusion coefficient (in _[L<sup>2</sup>/T]_) or gradient moments (in _[L<sup>-1</sup>]_).

`Quantity` objects contain their value in the base SI unit and may be converted to a compatible unit. The following code sample summarizes these features.

```cpp
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    // Unit suffix
    auto const repetition_time = 500_ms;
    // Unit as an object
    auto const echo_time = 100*ms;
    
    // Suffix and objects, along with type inference
    auto const diffusion = 0.89_um*um/ms;
    
    // Value in SI base unit
    auto const length = 180_cm;
    double const length_in_meters = length.magnitude; // equals to 1.8
    
    // Value in a compatible unit
    auto const duration = 1_h;
    double const duration_in_seconds = duration.convert_to(s); // equals to 3600
}
```
