#ifndef _faa6a046_30f6_4e87_91a6_033e2330b405
#define _faa6a046_30f6_4e87_91a6_033e2330b405

#include <array>
#include <tuple>
#include <utility>
#include <vector>

#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

namespace epg
{

namespace operators
{

/**
 * @brief Return the row-wise matrix corresponding to the single-pool EPG pulse
 * operator.
 */
SYCOMORE_API std::vector<Complex> pulse_single_pool(Real angle, Real phase);

/**
 * @brief Return the row-wise matrix corresponding to the two-pools exchange EPG
 * pulse operator.
 */
SYCOMORE_API
std::vector<Complex>
pulse_exchange(Real angle_a, Real phase_a, Real angle_b, Real phase_b);

/**
 * @brief Return the row-wise matrix corresponding to the two-pools 
 * magnetization transfer EPG pulse operator.
 */
SYCOMORE_API
std::vector<Complex>
pulse_magnetization_transfer(
    Real angle_a, Real phase_a, Real saturation_rate, Real duration);

/**
 * @brief Return the scalars associated respectively with relaxation of 
 * the F̃ states and Z̃ states.
 */
SYCOMORE_API std::pair<Real, Real> relaxation_single_pool(
    Real R1, Real R2, Real duration);

/**
 * @brief Return the exchange-relaxation matrices.
 */
SYCOMORE_API
std::pair<std::array<Complex, 8>, std::array<Real, 4>>
relaxation_exchange(
    Real R1_a, Real R2_a, Real R1_b, Real R2_b,
    Real k_a, Real k_b, Real delta_b,
    Real M0_a, Real M0_b,
    Real duration);

SYCOMORE_API
std::array<Complex, 8>
relaxation_and_exchange_T(
    Real R2_a, Real R2_b, Real k_a, Real k_b, Real delta_b);

SYCOMORE_API
std::array<Real, 4>
relaxation_and_exchange_L(Real R1_a, Real R1_b, Real k_a, Real k_b);

/**
 * @brief Return the scalars associated respectively with diffusion of 
 * respectively the F̃_k, F̃^*_{-k}, and Z̃_k states.
 */
template<typename T>
std::tuple<T, T, T> diffusion(Real D, Real duration, T const & k, Real delta_k);

/**
 * @brief Return the rotation expressed as a complex exponential associated
 * with phase accumulation of respectively the F̃_k and F̃^*_{-k} states.
 */
SYCOMORE_API std::pair<Complex, Complex> phase_accumulation(Real angle);

/**
 * @brief Return the phase accumulation expressed as a complex exponential 
 * associated with bulk motion of respectively the F̃_k, F̃^*_{-k}, and Z̃_k 
 * states.
 *
 * @param v: projection of the velocity vector on the axis associated with k.
 */
template<typename TReal, typename TComplex>
std::tuple<TComplex, TComplex, TComplex>
bulk_motion(Real v, Real duration, TReal const & k, Real delta_k);

}
    
}

}

#include "operators.txx"

#endif // _faa6a046_30f6_4e87_91a6_033e2330b405
