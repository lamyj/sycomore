#ifndef _d3b9de75_d62d_4445_b032_a08F985a5d10
#define _d3b9de75_d62d_4445_b032_a08F985a5d10

#include <utility>
#include <vector>

#include "sycomore/simd.h"
#include "sycomore/sycomore.h"

namespace sycomore 
{

namespace epg
{

namespace simd_api
{

#define SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(return_, name, parameters) \
    template<int InstructionSet> return_ name parameters; \
    template<> return_ name<0> parameters; \
    extern template return_ name<XSIMD_X86_SSE2_VERSION> parameters; \
    extern template return_ name<XSIMD_X86_AVX_VERSION> parameters; \
    extern template return_ name<XSIMD_X86_AVX512_VERSION> parameters;

// Functions with a _w suffix are worker functions, functions with a _d suffix
// are dispatcher functions.

/*******************************************************************************
 *                                Pulse operator                               *
 ******************************************************************************/

template<typename ValueType>
void apply_pulse_w(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, apply_pulse_d, 
    (
        std::vector<Complex> const & T, 
        Complex * F, Complex * F_star, Complex * Z,
        unsigned int states_count))

/*******************************************************************************
 *                             Relaxation operator                             *
 ******************************************************************************/

template<typename ValueType>
void relaxation_w(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, relaxation_d, 
    (
        std::pair<Real, Real> const & E,
        Real * F, Real * F_star, Real * Z, 
        unsigned int states_count))

/*******************************************************************************
 *                             Diffusion operator                              *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void diffusion_w(
    Real delta_k, Real tau, Real D, Real const * k_array,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, diffusion_d, 
    (
        Real delta_k, Real tau, Real D, Real const * k_array,
        Complex * F, Complex * F_star, Complex * Z,
        unsigned int states_count))

/*******************************************************************************
 *                           3D diffusion operator                             *
 ******************************************************************************/

template<typename ValueType>
void diffusion_3d_b_w(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D, 
    std::size_t begin, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, diffusion_3d_b_d, 
    (
        Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
        Real delta_k_product_term, Real tau, Real D_mn,
        Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D, 
        unsigned int states_count))

template<typename RealType, typename ComplexType>
void diffusion_3d_w(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, diffusion_3d_d, 
    (
        Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
        Complex * F, Complex * F_star, Complex * Z,
        unsigned int states_count))

/*******************************************************************************
 *                           Off-resonance operator                            *
 ******************************************************************************/

template<typename ValueType>
void off_resonance_w(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, off_resonance_d, 
    (
        std::pair<Complex, Complex> const & phi,
        Complex * F, Complex * F_star, Complex * Z, 
        unsigned int states_count))

/*******************************************************************************
 *                            Bulk motion operator                             *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void bulk_motion_w(
    Real delta_k, Real v, Real tau, Real const * k_array,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, bulk_motion_d, 
    (
        Real delta_k, Real v, Real tau, Real const * k_array,
        Complex * F, Complex * F_star, Complex * Z,
        unsigned int states_count))

/*******************************************************************************
 *                          Function table and set-up                          *
 ******************************************************************************/

extern decltype(&apply_pulse_d<0>) apply_pulse;
extern decltype(&relaxation_d<0>) relaxation;
extern decltype(&diffusion_d<0>) diffusion;
extern decltype(&diffusion_3d_b_d<0>) diffusion_3d_b;
extern decltype(&diffusion_3d_d<0>) diffusion_3d;
extern decltype(&off_resonance_d<0>) off_resonance;
extern decltype(&bulk_motion_d<0>) bulk_motion;

void set_api(int instruction_set);

bool set_default_api();

extern bool const api_is_set;

}

}

}

#include "simd_api.txx"

#endif // _d3b9de75_d62d_4445_b032_a08F985a5d10
