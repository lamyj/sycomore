#include <iostream>

#include <sycomore/epg/Discrete3D.h>
#include <sycomore/Species.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

#include <xtensor/xbuilder.hpp>

using namespace sycomore::units;

sycomore::Real dw_se(
    sycomore::Species const & species, sycomore::Quantity const & TE,
    sycomore::Vector3Q const & G)
{
    sycomore::epg::Discrete3D model(species);
    
    // Excitation
    model.apply_pulse(90*deg);
    
    // Diffusion-sensitization gradient is on all the time before the refocussing
    model.apply_time_interval(TE/2, G);
    
    // Refocussing
    model.apply_pulse(180*deg);
    
    // Diffusion-sensitization gradient is on all the time before TE
    model.apply_time_interval(TE/2, G);
    
    return std::abs(model.echo());
}

int main()
{
    sycomore::Species species(
        1000*ms, 100*ms,
        xt::diag(sycomore::Vector3R{2000, 500, 100})*std::pow(um, 2)/s);
    
    auto TE = 50*ms, b = 2000*s/std::pow(mm, 2);
    
    // Assuming Δ = δ (i.e. we can neglect the duration of the RF pulse, the
    // spatial encoding gradients and the crusher gradients, and the diffusion
    // gradient is on all the time, we get b = (γGδ)² ⋅ (2δ/3) with δ=TE/2.
    auto G = std::pow(
        b/(std::pow(sycomore::gamma, 2) * 2/3 * std::pow(TE/2, 3)), 0.5);
    std::cout << "G=" << G.convert_to(mT/m) << " mT/m\n";
    
    // Simulate the sequence with a diffusion gradient along different directions
    auto S_0 = dw_se(species, TE, {0*mT/m, 0*mT/m, 0*mT/m});
    auto S_1 = dw_se(species, TE, G*sycomore::Vector3R{1, 0, 0});
    auto S_2 = dw_se(species, TE, G*sycomore::Vector3R{0, 1, 0});
    auto S_3 = dw_se(species, TE, G*sycomore::Vector3R{0, 0, 1});
    
    std::cout
        << "S_x=" << S_1/S_0 << "\n"
        << "S_y=" << S_2/S_0 << "\n"
        << "S_z=" << S_3/S_0 << "\n";
    
    return 0;
}
