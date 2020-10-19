#ifndef _d3b9de75_d62d_4445_b032_a08F985a5d10
#define _d3b9de75_d62d_4445_b032_a08F985a5d10

#include <utility>
#include <vector>

#include "sycomore/sycomore.h"

namespace sycomore 
{

namespace epg
{

namespace simd_api
{

/*******************************************************************************
 *                                Pulse operator                               *
 ******************************************************************************/

template<typename ValueType>
void apply_pulse_d(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t start, std::size_t end, std::size_t step);

template<int InstructionSet>
void
apply_pulse_d(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

template<>
void
apply_pulse_d<0>(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

/*******************************************************************************
 *                             Relaxation operator                             *
 ******************************************************************************/

template<typename ValueType>
void relaxation_d(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z,
    std::size_t start, std::size_t end, std::size_t step);

template<int InstructionSet>
void
relaxation_d(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z, unsigned int states_count);

template<>
void
relaxation_d<0>(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z, unsigned int states_count);

/*******************************************************************************
 *                             Diffusion operator                              *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void diffusion_d(
    Real delta_k, Real tau, Real D, Real const * k_array,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

template<int InstructionSet>
void
diffusion_d(
    Real delta_k, Real tau, Real D, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

template<>
void
diffusion_d<0>(
    Real delta_k, Real tau, Real D, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

/*******************************************************************************
 *                           3D diffusion operator                             *
 ******************************************************************************/

template<typename ValueType>
void diffusion_3d_b_d(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D, 
    std::size_t begin, std::size_t end, std::size_t step);

template<int InstructionSet>
void
diffusion_3d_b_d(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D,
    unsigned int states_count);

template<>
void
diffusion_3d_b_d<0>(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D,
    unsigned int states_count);

template<typename RealType, typename ComplexType>
void diffusion_3d_d(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

template<int InstructionSet>
void
diffusion_3d_d(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

template<>
void
diffusion_3d_d<0>(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

/*******************************************************************************
 *                           Off-resonance operator                            *
 ******************************************************************************/

template<typename ValueType>
void off_resonance_d(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

template<int InstructionSet>
void
off_resonance_d(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

template<>
void
off_resonance_d<0>(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

/*******************************************************************************
 *                            Bulk motion operator                             *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void bulk_motion_d(
    Real delta_k, Real v, Real tau, Real const * k_array,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step);

template<int InstructionSet>
void
bulk_motion_d(
    Real delta_k, Real v, Real tau, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

template<>
void
bulk_motion_d<0>(
    Real delta_k, Real v, Real tau, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

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
