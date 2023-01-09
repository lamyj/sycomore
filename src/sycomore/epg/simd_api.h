#ifndef _d3b9de75_d62d_4445_b032_a08F985a5d10
#define _d3b9de75_d62d_4445_b032_a08F985a5d10

#include <utility>
#include <vector>

#include "sycomore/epg/Model.h"
#include "sycomore/simd.h"
#include "sycomore/sycomore.h"

namespace sycomore 
{

namespace epg
{

namespace simd_api
{

#define SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(return_, name, parameters) \
    template<INSTRUCTION_SET_TYPE InstructionSet> return_ name parameters; \
    template<> return_ name<unsupported> parameters; \
    extern template return_ name<XSIMD_X86_SSE2_VERSION> parameters; \
    extern template return_ name<XSIMD_X86_AVX_VERSION> parameters; \
    extern template return_ name<XSIMD_X86_AVX512_VERSION> parameters;

// Functions with a _w suffix are worker functions, functions with a _d suffix
// are dispatcher functions.

/*******************************************************************************
 *                                Pulse operators                              *
 ******************************************************************************/

template<typename ValueType>
void apply_pulse_single_pool_w(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
    void, apply_pulse_single_pool_d, 
    (std::vector<Complex> const & T, Model & model, std::size_t states_count))

template<typename ValueType>
void apply_pulse_exchange_w(
    std::vector<Complex> const & T,
    Complex * F_a, Complex * F_star_a, Complex * Z_a,
    Complex * F_b, Complex * F_star_b, Complex * Z_b,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
    void, apply_pulse_exchange_d, 
    (std::vector<Complex> const & T, Model & model, std::size_t states_count))

template<typename ValueType>
void apply_pulse_magnetization_transfer_w(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z_a, Complex * Z_b,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
    void, apply_pulse_magnetization_transfer_d, 
    (std::vector<Complex> const & T, Model & model, std::size_t states_count))

/*******************************************************************************
 *                             Relaxation operators                            *
 ******************************************************************************/

template<typename ValueType>
void relaxation_single_pool_w(
    std::pair<Real, Real> const & E,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
    void, relaxation_single_pool_d, 
    (std::pair<Real, Real> const & E, Model & model, std::size_t states_count))

template<typename ValueType>
void relaxation_exchange_w(
    std::array<Complex, 8> const & Xi_T, std::array<Real, 4> const & Xi_L,
    Complex * F_a, Complex * F_star_a, Complex * Z_a,
    Complex * F_b, Complex * F_star_b, Complex * Z_b,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, relaxation_exchange_d, 
    (
        std::array<Complex, 8> const & Xi_T, std::array<Real, 4> const & Xi_L,
        Model & model, std::size_t states_count))

template<typename ValueType>
void relaxation_magnetization_transfer_w(
    Real const & Xi_T, std::array<Real, 4> const & Xi_L,
    Complex * F_a, Complex * F_star_a, Complex * Z_a,
    Complex * F_b, Complex * F_star_b, Complex * Z_b,
    std::size_t start, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
   void, relaxation_magnetization_transfer_d, 
    (
        Real const & Xi_T, std::array<Real, 4> const & Xi_L,
        Model & model, std::size_t states_count))

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
        Model::Population & F, Model::Population & F_star, Model::Population & Z,
        std::size_t states_count))

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
    Complex * F, Complex * F_star,
    std::size_t begin, std::size_t end, std::size_t step);

SYCOMORE_DEFINE_SIMD_DISPATCHER_FUNCTION(
    void, off_resonance_d, 
    (
        std::pair<Complex, Complex> const & phi, 
        Model::Population & F, Model::Population & F_star,
        std::size_t states_count))

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
        Real delta_k, Real v, Real tau, Real const * k_array, Model & model,
        std::size_t states_count))

/*******************************************************************************
 *                          Function table and set-up                          *
 ******************************************************************************/

extern decltype(&apply_pulse_single_pool_d<unsupported>) apply_pulse_single_pool;
extern decltype(&apply_pulse_exchange_d<unsupported>) apply_pulse_exchange;
extern decltype(&apply_pulse_magnetization_transfer_d<unsupported>)
    apply_pulse_magnetization_transfer;
extern decltype(&relaxation_single_pool_d<unsupported>) relaxation_single_pool;
extern decltype(&relaxation_exchange_d<unsupported>) relaxation_exchange;
extern decltype(&relaxation_magnetization_transfer_d<unsupported>)
    relaxation_magnetization_transfer;
extern decltype(&diffusion_d<unsupported>) diffusion;
extern decltype(&diffusion_3d_b_d<unsupported>) diffusion_3d_b;
extern decltype(&diffusion_3d_d<unsupported>) diffusion_3d;
extern decltype(&off_resonance_d<unsupported>) off_resonance;
extern decltype(&bulk_motion_d<unsupported>) bulk_motion;

void set_api(unsigned instruction_set);

bool set_default_api();

extern bool const api_is_set;

}

}

}

#include "simd_api.txx"

#endif // _d3b9de75_d62d_4445_b032_a08F985a5d10
