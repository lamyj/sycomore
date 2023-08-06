#include <iostream>

#include <sycomore/Species.h>
#include <sycomore/units.h>

#include <xtensor/xio.hpp>

int main()
{
    using namespace sycomore::units;
    
    sycomore::Species species(1000*ms, 100*ms);
    // Assign the diffusion coefficient as a scalar
    species.set_D(3*std::pow(um, 2)/s);
    // The diffusion coefficient is stored on the diagonal of the tensor
    std::cout << species.D()(0,0) << "\n";
    
    // Assign the diffusion coefficient as a tensor
    species.set_D( {
         {3*std::pow(um, 2)/s, 0*std::pow(um, 2)/s, 0*std::pow(um, 2)/s },
         {0*std::pow(um, 2)/s, 2*std::pow(um, 2)/s, 0*std::pow(um, 2)/s },
         {0*std::pow(um, 2)/s, 0*std::pow(um, 2)/s, 1*std::pow(um, 2)/s } });
    std::cout << species.D() << "\n";
    
    return 0;
}
