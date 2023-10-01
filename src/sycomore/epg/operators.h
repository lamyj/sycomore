#ifndef _faa6a046_30f6_4e87_91a6_033e2330b405
#define _faa6a046_30f6_4e87_91a6_033e2330b405

#include <array>
#include <tuple>
#include <utility>
#include <vector>

#include "sycomore/sycomore.h"

namespace sycomore
{

namespace epg
{

namespace operators
{

/// @addtogroup EPGOperators
/// @{

/**
 * @brief Return the row-wise matrix corresponding to the single-pool EPG pulse
 * operator.
 */
std::vector<Complex> pulse_single_pool(Real angle, Real phase);

/**
 * @brief Return the row-wise matrix corresponding to the two-pools exchange EPG
 * pulse operator.
 */
std::vector<Complex>
pulse_exchange(Real angle_a, Real phase_a, Real angle_b, Real phase_b);

/**
 * @brief Return the row-wise matrix corresponding to the two-pools 
 * magnetization transfer EPG pulse operator.
 */
std::vector<Complex>
pulse_magnetization_transfer(Real angle_a, Real phase_a, Real saturation);

/**
 * @brief Return the scalars associated respectively with relaxation of 
 * the \f$\tilde{F}\f$ states and \f$\tilde{Z}\f$ states.
 */
std::pair<Real, Real> relaxation_single_pool(Real R1, Real R2, Real duration);

/**
 * @brief Return the exchange-relaxation matrices: non-zero terms of
 * \f$\Xi_T\f$ (row-major), row-major \f$\Xi_L\f$, and
 * \f$(\Xi_L - I) \Lambda_L^{-1} C\f$
 */
std::tuple<std::array<Complex, 8>, std::array<Real, 4>, std::array<Real, 2>>
relaxation_exchange(
    Real R1_a, Real R2_a, Real R1_b, Real R2_b,
    Real k_a, Real k_b, Real delta_b,
    Real M0_a, Real M0_b,
    Real duration);

/**
 * @brief Return the magnetization transfer relaxation matrices: 
 * \f$\Xi_T = e^{-R_{2a} \tau}\f$, row-major \f$\Xi_L\f$ and
 * \f$(\Xi_L - I) \Lambda_L^{-1} C\f$
 */
std::tuple<Real, std::array<Real, 4>, std::array<Real, 2>>
relaxation_magnetization_transfer(
    Real R1_a, Real R2_a, Real R1_b,
    Real k_a, Real k_b,
    Real M0_a, Real M0_b,
    Real duration);

/// @brief Return the non-zero terms of \f$\Xi_T\f$ (row-major)
std::array<Complex, 8>
relaxation_and_exchange_T(
    Real R2_a, Real R2_b, Real k_a, Real k_b, Real delta_b, Real duration);

/// @brief Return \f$\Xi_L\f$ (row-major)
std::array<Real, 4>
relaxation_and_exchange_L(
    Real R1_a, Real R1_b, Real k_a, Real k_b, Real duration);

/// @brief Return \f$(\Xi_L - I) \Lambda_L^{-1} C\f$
std::array<Real, 2>
relaxation_and_exchange_recovery(
    Real R1_a, Real R1_b, Real k_a, Real k_b, Real M0_a, Real M0_b,
    std::array<Real, 4> const & Xi_L);

/// @brief Closed-form 2x2 matrix exponential
std::array<Real, 4> expm(std::array<Real, 4> const & A);

/**
 * @brief Return the scalars associated respectively with diffusion of 
 * respectively the \f$\tilde{F}(k)\f$, \f$\tilde{F}^*(-k)\f$, and
 * \f$\tilde{Z}(k)\f$ states.
 */
template<typename T>
std::tuple<T, T, T> diffusion(Real D, Real duration, T const & k, Real delta_k);

/**
 * @brief Return the rotation expressed as a complex exponential associated
 * with phase accumulation of respectively the \f$\tilde{F}(k)\f$ and
 * \f$\tilde{F}^*(-k)\f$ states.
 */
std::pair<Complex, Complex> phase_accumulation(Real angle);

/**
 * @brief Return the phase accumulation expressed as a complex exponential 
 * associated with bulk motion of respectively the \f$\tilde{F}(k)\f$,
 * \f$\tilde{F}^*(-k)\f$, and \f$\tilde{Z}(k)\f$ states.
 *
 * @param v: projection of the velocity vector on the axis associated with k.
 */
template<typename TReal, typename TComplex>
std::tuple<TComplex, TComplex, TComplex>
bulk_motion(Real v, Real duration, TReal const & k, Real delta_k);

/// @}

}
    
}

}

#include "operators.txx"

#endif // _faa6a046_30f6_4e87_91a6_033e2330b405
