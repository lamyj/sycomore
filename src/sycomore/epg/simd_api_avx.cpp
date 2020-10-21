#include "simd_api.h"

namespace sycomore
{

namespace epg
{

namespace simd_api
{

template 
void apply_pulse_d<XSIMD_X86_AVX_VERSION>(
    std::vector<Complex> const & T, 
    Complex * F, Complex * F_star, Complex * Z,
    unsigned int states_count);

template 
void relaxation_d<XSIMD_X86_AVX_VERSION>(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z, unsigned int states_count);

template 
void diffusion_d<XSIMD_X86_AVX_VERSION>(
    Real delta_k, Real tau, Real D, Real const * k_array,
    Complex * F, Complex * F_star, Complex * Z,
    unsigned int states_count);

template 
void diffusion_3d_b_d<XSIMD_X86_AVX_VERSION>(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D, 
    unsigned int states_count);

template
void diffusion_3d_d<XSIMD_X86_AVX_VERSION>(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z,
    unsigned int states_count);

template
void
off_resonance_d<XSIMD_X86_AVX_VERSION>(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

template
void
bulk_motion_d<XSIMD_X86_AVX_VERSION>(
    Real delta_k, Real v, Real tau, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count);

}

}

}
