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
    using std::cos; using std::exp; using std::pow;
    
    auto const a = angle.convert_to(units::rad);
    auto const p = phase.convert_to(units::rad);
    
    constexpr Complex const i{0,1};

    return std::vector<Complex>{
        pow(cos(a/2), 2),              exp(2.*i*p)*pow(sin(a/2), 2), -i*exp( i*p)*sin(a),
        exp(-2.*i*p)*pow(sin(a/2), 2), pow(cos(a/2), 2),              i*exp(-i*p)*sin(a),
        -i/2.*exp(-i*p)*sin(a),        i/2.*exp(i*p)*sin(a),          cos(a)
    };
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
    using std::exp; using std::pow;
    
    auto const k = k_.magnitude;
    auto const tau = duration.magnitude;
    auto const delta_k = delta_k_.magnitude;
    auto const d = species.get_D()[0].magnitude;
        
    auto const b_T_plus = tau*(pow(k+delta_k/2, 2) + pow(delta_k, 2) / 12);
    auto const D_T_plus = exp(-b_T_plus*d);
    
    auto const b_T_minus = tau*(pow(-k+delta_k/2, 2) + pow(delta_k, 2) / 12);
    auto const D_T_minus = exp(-b_T_minus*d);
    
    auto const b_L = pow(k, 2) * tau;
    auto const D_L = exp(-b_L*d);
    
    return std::make_tuple(D_T_plus, D_T_minus, D_L);
}

}
    
}

}
