#include "operators.h"

#include <tuple>
#include <utility>
#include <vector>

#include "sycomore/Quantity.h"
#include "sycomore/Species.h"

namespace sycomore
{

namespace epg
{

namespace operators
{

std::vector<Complex> pulse(Quantity angle, Quantity phase)
{
    auto const a = angle.convert_to(units::rad);
    auto const p = phase.convert_to(units::rad);
    
    constexpr Complex const i{0,1};
    
    std::vector<Complex> T(9, 0);
    T[3*0+0] = std::pow(std::cos(a/2), 2);
    T[3*0+1] = std::exp(2.*i * p)*std::pow(std::sin(a/2), 2);
    T[3*0+2] = -i*std::exp(i*p)*std::sin(a);
    T[3*1+0] = std::exp(-2.*i*p)*std::pow(std::sin(a/2), 2);
    T[3*1+1] = std::pow(std::cos(a/2), 2);
    T[3*1+2] = i*std::exp(-i*p)*std::sin(a);
    T[3*2+0] = -i/2.*std::exp(-i*p)*std::sin(a);
    T[3*2+1] = i/2.*std::exp(i*p)*std::sin(a);
    T[3*2+2] = std::cos(a);
    
    return T;
}

std::pair<Real, Real> 
relaxation(Species const & species, Quantity const & duration)
{
    auto const E_1 = std::exp((-duration*species.get_R1()).magnitude);
    auto const E_2 = std::exp((-duration*species.get_R2()).magnitude);
    return std::make_pair(E_1, E_2);
}

std::tuple<Real, Real, Real>
diffusion(
    Species const & species, Quantity const & duration, 
    Quantity const & k_, Quantity const & delta_k_)
{
    // NOTE: b_T differs between F̃^+ and F̃^{-*} since F̃^{-*}(k) is F(-k^*)
    
    using namespace units;
    
    auto const k = k_.magnitude;
    auto const tau = duration.magnitude;
    auto const delta_k = delta_k_.magnitude;
    auto const d = species.get_D().magnitude;
        
    auto const b_T_plus = tau*(std::pow(k+delta_k/2, 2) + std::pow(delta_k, 2) / 12);
    auto const D_T_plus = std::exp(-b_T_plus*d);
    
    auto const b_T_minus = tau*(std::pow(-k+delta_k/2, 2) + std::pow(delta_k, 2) / 12);
    auto const D_T_minus = std::exp(-b_T_minus*d);
    
    auto const b_L = std::pow(k, 2) * tau;
    auto const D_L = std::exp(-b_L*d);
    
    return std::make_tuple(D_T_plus, D_T_minus, D_L);
}

}
    
}

}
