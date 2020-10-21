#ifndef _9d63e9db_c0d8_48c2_b55a_549aa72a0440
#define _9d63e9db_c0d8_48c2_b55a_549aa72a0440

#include "operators.h"

#include <cmath>
#include <tuple>
#include <utility>
#include <vector>

#include "sycomore/simd.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

namespace operators
{

template<typename T>
std::tuple<T, T, T> diffusion(
    Real D, Real duration, T const & k, Real delta_k)
{
    auto const b_T_plus = 
        duration*(sycomore::simd::pow(k+delta_k/2, 2.) + std::pow(delta_k, 2.) / 12);
    auto const D_T_plus = sycomore::simd::exp(-b_T_plus*D);
    
    auto const b_T_minus = 
        duration*(sycomore::simd::pow(-k+delta_k/2, 2.) + std::pow(delta_k, 2.) / 12);
    auto const D_T_minus = sycomore::simd::exp(-b_T_minus*D);
    
    auto const b_L = sycomore::simd::pow(k, 2.) * duration;
    auto const D_L = sycomore::simd::exp(-b_L*D);
    
    return std::make_tuple(D_T_plus, D_T_minus, D_L);
}

template<typename TReal, typename TComplex>
std::tuple<TComplex, TComplex, TComplex>
bulk_motion(Real v, Real duration, TReal const & k, Real delta_k)
{
    // NOTE: the following form matches both batch and non-batch complex.
    static TComplex const i{Complex{0., 1.}};
    
    auto const v_tau = v*duration;
    
    // FIXME? In Weigel 2015, eq. 48 gives the phase accumulation as k⋅v⋅Δt and
    // states
    // > This phase term has to be added to the phase of a given configuration 
    // > state by the coherent motion operator.
    // However, there is a sign inversion in the following steps
    // > This is best done by multiplying a complex phase factor of the type 
    // > exp(-iΔϕ)
    // The same thing appears in Sodickson 1998.
    auto const J_T = sycomore::simd::exp(i * (k+delta_k/2) * v_tau);
    auto const J_L = sycomore::simd::exp(i * k * v_tau);
    
    return std::make_tuple(J_T, sycomore::simd::conj(J_T), J_L);
}

}

}

}

#endif // _9d63e9db_c0d8_48c2_b55a_549aa72a0440
