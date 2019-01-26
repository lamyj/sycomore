# Species

A species is characterized by its relaxation rates (R<sub>1</sub>, R<sub>2</sub> and R<sub>2</sub>'), its diffusivity D and its relative resonance frequency Δ<sub>ω</sub>. It is described by the `Species` class, which stores the species parameters in their base SI units (i.e. _Hz_ for R<sub>1</sub>, R<sub>2</sub> and R<sub>2</sub>', _m<sup>2</sup>/s_ for D and _rad/s_ for Δ<sub>ω</sub>). R<sub>1</sub>, R<sub>2</sub> are mandatory parameters, all others are optional and are equal to 0 if unspecified. Using the [units system](units.md), `Species` objects can be created using different syntaxes. 

The parameters of the constructors of the `Species` class are (in order):

- `R1` (mandatory)
- `R2` (mandatory)
- `D` (optional, defaults to 0)
- `R2_prime` (optional, defaults to 0)
- `delta_omega` (optional, defaults to 0)

```cpp
#include <sycomore/Species.h>
#include <sycomore/units.h>

void main()
{
    // In this example, the relaxation paramters are always the same, but 
    // specified using different syntaxes.
    
    using namespace sycomore::units;
    
    // Values in base SI values, default parameters (D=0, R2_prime=0, delta_omega=0)
    sycomore::Species const foo(1.0, 0.1);
    
    // Values using times, user-specified diffusivity
    sycomore::Species const bar(1000_ms, 100_ms, 0.89_um*um/s);
    
    // Values using rates
    sycomore::Species const bar(1_Hz, 1/100_ms);
}
```
