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

std::vector<Complex> pulse(Real angle, Real phase)
{
    using std::cos; using std::exp; using std::pow;
    
    auto const & a = angle;
    auto const & p = phase;
    
    constexpr Complex const i{0,1};

    return {
        pow(cos(a/2), 2),              exp(2.*i*p)*pow(sin(a/2), 2), -i*exp( i*p)*sin(a),
        exp(-2.*i*p)*pow(sin(a/2), 2), pow(cos(a/2), 2),              i*exp(-i*p)*sin(a),
        -i/2.*exp(-i*p)*sin(a),        i/2.*exp(i*p)*sin(a),          cos(a)
    };
}

std::pair<Real, Real> 
relaxation(Real R1, Real R2, Real duration)
{
    auto const E_1 = std::exp(-duration*R1);
    auto const E_2 = std::exp(-duration*R2);
    return std::make_pair(E_1, E_2);
}

std::tuple<Real, Real, Real>
diffusion(Real D, Real duration, Real k, Real delta_k)
{
    // NOTE: b_T differs between F̃^+ and F̃^{-*} since F̃^{-*}(k) is F(-k^*)
    
    using std::exp; using std::pow;
    
    auto const b_T_plus = duration*(pow(k+delta_k/2, 2.) + pow(delta_k, 2.) / 12);
    auto const D_T_plus = exp(-b_T_plus*D);
    
    auto const b_T_minus = duration*(pow(-k+delta_k/2, 2.) + pow(delta_k, 2.) / 12);
    auto const D_T_minus = exp(-b_T_minus*D);
    
    auto const b_L = pow(k, 2.) * duration;
    auto const D_L = exp(-b_L*D);
    
    return std::make_tuple(D_T_plus, D_T_minus, D_L);
}

std::pair<Complex, Complex> phase_accumulation(Real angle)
{
    constexpr Complex const i{0,1};
    return {std::exp(i*angle), std::exp(i*-angle)};
}

std::tuple<Complex, Complex, Complex>
bulk_motion(Real v, Real duration, Real k, Real delta_k)
{
    constexpr Complex const i{0,1};
    
    auto const v_tau = v*duration;
    
    // FIXME? In Weigel 2015, eq. 48 gives the phase accumulation as k⋅v⋅Δt and
    // states
    // > This phase term has to be added to the phase of a given configuration 
    // > state by the coherent motion operator.
    // However, there is a sign inversion in the following steps
    // > This is best done by multiplying a complex phase factor of the type 
    // > exp(-iΔϕ)
    // The same thing appears in Sodickson 1998.
    auto const J_T = std::exp(i*double((k+delta_k/2)*v_tau));
    auto const J_L = std::exp(i*double(k*v_tau));
    
    return std::make_tuple(J_T, std::conj(J_T), J_L);
}

}
    
}

}
