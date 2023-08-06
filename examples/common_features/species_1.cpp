#include <sycomore/Species.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    // Create a Species from either relaxation times, relaxation rates or both
    sycomore::Species species_1(1000*ms, 100*ms);
    sycomore::Species species_2(1*Hz, 10*Hz);
    sycomore::Species species_3(1000*ms, 10*Hz);
    
    return 0;
}
