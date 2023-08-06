#include <sycomore/epg/Discrete.h>
#include <sycomore/Species.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    sycomore::Species species_a(779*ms, 45*ms);
    sycomore::epg::Discrete single_pool_model(species_a);
    
    sycomore::Species species_b(100*ms, 20*ms);
    sycomore::Vector3R M0_a{0, 0, 0.8}, M0_b{0, 0, 0.2};
    auto k_a = 2*Hz;
    sycomore::epg::Discrete exchange_model(species_a, species_b, M0_a, M0_b, k_a);
    
    auto R1_b = 100*ms;
    sycomore::epg::Discrete mt_model(species_a, R1_b, M0_a, M0_b, k_a);
}
