#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    // Combination of base units
    auto diffusion_coefficient = 0.89*std::pow(um, 2)/ms;
    
    // Magnitude of the quantity, in SI unit
    auto length = 180*cm;
    auto length_in_meters = length.magnitude; // Equals to 1.8
    
    // Conversion: magnitude of the quantity, in specified unit
    auto duration = 1*h;
    auto duration_in_seconds = duration.convert_to(s); // Equals to 3600
    
    return 0;
}
