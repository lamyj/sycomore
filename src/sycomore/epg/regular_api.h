#ifndef _d3b9de75_d62d_4445_b032_a08F985a5d10
#define _d3b9de75_d62d_4445_b032_a08F985a5d10

#include <utility>
#include <vector>

#include "sycomore/simd.h"

namespace sycomore 
{

namespace epg
{

namespace regular_api
{

/*******************************************************************************
 *                                Pulse operator                               *
 ******************************************************************************/

template<typename ValueType>
void apply_pulse_d(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t start, std::size_t end, std::size_t step)
{
    for(std::size_t i=start; i<end; i+=step)
    {
        ValueType F_i, F_star_i, Z_i;
        sycomore::simd::load_aligned(F+i, F_i);
        sycomore::simd::load_aligned(F_star+i, F_star_i);
        sycomore::simd::load_aligned(Z+i, Z_i);
    
        auto const F_i_new = 
            T[3*0+0] * F_i + T[3*0+1] * F_star_i + T[3*0+2] * Z_i;
        auto const F_i_star_new = 
            T[3*1+0] * F_i + T[3*1+1] * F_star_i + T[3*1+2] * Z_i;
        auto const Z_i_new =
            T[3*2+0] * F_i + T[3*2+1] * F_star_i + T[3*2+2] * Z_i;
        
        sycomore::simd::store_aligned(F_i_new, F+i);
        sycomore::simd::store_aligned(F_i_star_new, F_star+i);
        sycomore::simd::store_aligned(Z_i_new, Z+i);
    }
}

template<int InstructionSet>
void
apply_pulse_d(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using Batch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % Batch::size;
    
    apply_pulse_d<Batch>(T, F, F_star, Z, 0, simd_end, Batch::size);
    apply_pulse_d<Complex>(T, F, F_star, Z, simd_end, states_count, 1);
}

template<>
void
apply_pulse_d<0>(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    apply_pulse_d<Complex>(T, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                             Relaxation operator                             *
 ******************************************************************************/

template<typename ValueType>
void relaxation_d(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z,
    std::size_t start, std::size_t end, std::size_t step)
{
    for(std::size_t i=start; i<end; i+=step)
    {
        ValueType F_i, F_star_i, Z_i;
        sycomore::simd::load_aligned(F+i, F_i);
        sycomore::simd::load_aligned(F_star+i, F_star_i);
        sycomore::simd::load_aligned(Z+i, Z_i);
        
        sycomore::simd::store_aligned(F_i*E.second, F+i);
        sycomore::simd::store_aligned(F_star_i*E.second, F_star+i);
        sycomore::simd::store_aligned(Z_i*E.first, Z+i);
    }
}

template<int InstructionSet>
void
relaxation_d(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z, unsigned int states_count)
{
    using Batch = simd::Batch<Real, InstructionSet>;
    // Use 2*states_count as we are getting reinterpreted real arrays
    auto const simd_end = 2*states_count - 2*states_count % Batch::size;
    
    relaxation_d<Batch>(E, F, F_star, Z, 0, simd_end, Batch::size);
    relaxation_d<Real>(E, F, F_star, Z, simd_end, 2*states_count, 1);
}

template<>
void
relaxation_d<0>(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z, unsigned int states_count)
{
    // Pass 2*states_count as we are getting reinterpreted real arrays
    relaxation_d<Real>(E, F, F_star, Z, 0, 2*states_count, 1);
}

/*******************************************************************************
 *                             Diffusion operator                              *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void diffusion_d(
    Real delta_k, Real tau, Real D, Real const * k_array,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step)
{
    for(std::size_t i=begin; i<end; i+=step)
    {
        RealType k;
        sycomore::simd::load_aligned(k_array+i, k);
        
        auto const D_ = operators::diffusion(D, tau, k, delta_k);
        
        ComplexType F_i;
        sycomore::simd::load_aligned(F+i, F_i);
        sycomore::simd::store_aligned(F_i*std::get<0>(D_), F+i);
        
        ComplexType F_star_i;
        sycomore::simd::load_aligned(F_star+i, F_star_i);
        sycomore::simd::store_aligned(F_star_i*std::get<1>(D_), F_star+i);
        
        ComplexType Z_i;
        sycomore::simd::load_aligned(Z+i, Z_i);
        sycomore::simd::store_aligned(Z_i*std::get<2>(D_), Z+i);
    }
}

template<int InstructionSet>
void
diffusion_d(
    Real delta_k, Real tau, Real D, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using RealBatch = simd::Batch<Real, InstructionSet>;
    using ComplexBatch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % ComplexBatch::size;
    
    diffusion_d<RealBatch, ComplexBatch>(
        delta_k, tau, D, k, F, F_star, Z, 0, simd_end, ComplexBatch::size);
    diffusion_d<Real, Complex>(
        delta_k, tau, D, k, F, F_star, Z, simd_end, states_count, 1);
}

template<>
void
diffusion_d<0>(
    Real delta_k, Real tau, Real D, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    diffusion_d<Real, Complex>(
        delta_k, tau, D, k, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                           Off-resonance operator                            *
 ******************************************************************************/

template<typename ValueType>
void off_resonance_d(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step)
{
    for(int i=begin; i<end; i+=step)
    {
        ValueType F_i, F_star_i;
        sycomore::simd::load_aligned(F+i, F_i);
        sycomore::simd::load_aligned(F_star+i, F_star_i);
        
        sycomore::simd::store_aligned(F_i*phi.first, F+i);
        sycomore::simd::store_aligned(F_star_i*phi.second, F_star+i);
        
        // Z̃ states are unaffected
    }
}

template<int InstructionSet>
void
off_resonance_d(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using Batch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % Batch::size;
    
    off_resonance_d<Batch>(phi, F, F_star, Z, 0, simd_end, Batch::size);
    off_resonance_d<Complex>(phi, F, F_star, Z, simd_end, states_count, 1);
}

template<>
void
off_resonance_d<0>(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    off_resonance_d<Complex>(phi, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                            Bulk motion operator                             *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void bulk_motion_d(
    Real delta_k, Real v, Real tau, Real const * k_array,
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step)
{
    for(std::size_t i=begin; i<end; i+=step)
    {
        RealType k;
        sycomore::simd::load_aligned(k_array+i, k);
        
        auto const J = operators::bulk_motion<RealType, ComplexType>(
            v, tau, k, delta_k);
        
        ComplexType F_i;
        sycomore::simd::load_aligned(F+i, F_i);
        sycomore::simd::store_aligned(F_i*std::get<0>(J), F+i);
        
        ComplexType F_star_i;
        sycomore::simd::load_aligned(F_star+i, F_star_i);
        sycomore::simd::store_aligned(F_star_i*std::get<1>(J), F_star+i);
        
        ComplexType Z_i;
        sycomore::simd::load_aligned(Z+i, Z_i);
        sycomore::simd::store_aligned(Z_i*std::get<2>(J), Z+i);
    }
}

template<int InstructionSet>
void
bulk_motion_d(
    Real delta_k, Real v, Real tau, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using RealBatch = simd::Batch<Real, InstructionSet>;
    using ComplexBatch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % ComplexBatch::size;
    
    diffusion_d<RealBatch, ComplexBatch>(
        delta_k, v, tau, k, F, F_star, Z, 0, simd_end, ComplexBatch::size);
    diffusion_d<Real, Complex>(
        delta_k, v, tau, k, F, F_star, Z, simd_end, states_count, 1);
}

template<>
void
bulk_motion_d<0>(
    Real delta_k, Real v, Real tau, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    diffusion_d<Real, Complex>(
        delta_k, v, tau, k, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                          Function table and set-up                          *
 ******************************************************************************/

decltype(&apply_pulse_d<0>) apply_pulse = nullptr;
decltype(&relaxation_d<0>) relaxation = nullptr;
decltype(&diffusion_d<0>) diffusion = nullptr;
decltype(&off_resonance_d<0>) off_resonance = nullptr;
decltype(&bulk_motion_d<0>) bulk_motion = nullptr;

void set_api(int instruction_set)
{
    SYCOMORE_SET_API_FUNCTION(apply_pulse)
    SYCOMORE_SET_API_FUNCTION(relaxation)
    SYCOMORE_SET_API_FUNCTION(diffusion)
    SYCOMORE_SET_API_FUNCTION(off_resonance)
    SYCOMORE_SET_API_FUNCTION(bulk_motion)
}

bool set_default_api()
{
    set_api(simd::get_instruction_set());
    return true;
}

bool const api_is_set=set_default_api();

}

}

}

#endif // _d3b9de75_d62d_4445_b032_a08F985a5d10
