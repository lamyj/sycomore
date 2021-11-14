#ifndef _c90ddae6_b878_4da1_b3a4_36c025c09ae1
#define _c90ddae6_b878_4da1_b3a4_36c025c09ae1

#include "simd_api.h"

#include <utility>
#include <vector>

#include "sycomore/epg/operators.h"
#include "sycomore/simd.h"
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
void apply_pulse_w(
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

template<INSTRUCTION_SET_TYPE InstructionSet>
void
apply_pulse_d(
    std::vector<Complex> const & T,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using Batch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % Batch::size;
    
    apply_pulse_w<Batch>(T, F, F_star, Z, 0, simd_end, Batch::size);
    apply_pulse_w<Complex>(T, F, F_star, Z, simd_end, states_count, 1);
}

/*******************************************************************************
 *                             Relaxation operator                             *
 ******************************************************************************/

template<typename ValueType>
void relaxation_w(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z,
    std::size_t start, std::size_t end, std::size_t step)
{
    for(std::size_t i=start; i<end; i+=step)
    {
        ValueType F_i;
        sycomore::simd::load_aligned(F+i, F_i);
        sycomore::simd::store_aligned(F_i*E.second, F+i);

        ValueType F_star_i;
        sycomore::simd::load_aligned(F_star+i, F_star_i);
        sycomore::simd::store_aligned(F_star_i*E.second, F_star+i);

        ValueType Z_i;
        sycomore::simd::load_aligned(Z+i, Z_i);
        sycomore::simd::store_aligned(Z_i*E.first, Z+i);
    }
}

template<INSTRUCTION_SET_TYPE InstructionSet>
void
relaxation_d(
    std::pair<Real, Real> const & E,
    Real * F, Real * F_star, Real * Z, unsigned int states_count)
{
    using Batch = simd::Batch<Real, InstructionSet>;
    // Use 2*states_count as we are getting reinterpreted real arrays
    auto const simd_end = 2*states_count - 2*states_count % Batch::size;
    
    relaxation_w<Batch>(E, F, F_star, Z, 0, simd_end, Batch::size);
    relaxation_w<Real>(E, F, F_star, Z, simd_end, 2*states_count, 1);
}

/*******************************************************************************
 *                             Diffusion operator                              *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void diffusion_w(
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

template<INSTRUCTION_SET_TYPE InstructionSet>
void
diffusion_d(
    Real delta_k, Real tau, Real D, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using RealBatch = simd::Batch<Real, InstructionSet>;
    using ComplexBatch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % ComplexBatch::size;
    
    diffusion_w<RealBatch, ComplexBatch>(
        delta_k, tau, D, k, F, F_star, Z, 0, simd_end, ComplexBatch::size);
    diffusion_w<Real, Complex>(
        delta_k, tau, D, k, F, F_star, Z, simd_end, states_count, 1);
}

/*******************************************************************************
 *                           3D diffusion operator                             *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void diffusion_3d_w(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z,
    std::size_t begin, std::size_t end, std::size_t step)
{
    for(int i=begin; i<end; i+=step)
    {
        ComplexType F_i; RealType b_T_plus_D_i;
        sycomore::simd::load_aligned(F+i, F_i);
        sycomore::simd::load_aligned(b_T_plus_D+i, b_T_plus_D_i);
        sycomore::simd::store_aligned(F_i * simd::exp(-b_T_plus_D_i), F+i);
        
        ComplexType F_star_i; RealType b_T_minus_D_i;
        sycomore::simd::load_aligned(F_star+i, F_star_i);
        sycomore::simd::load_aligned(b_T_minus_D+i, b_T_minus_D_i);
        sycomore::simd::store_aligned(
            F_star_i * simd::exp(-b_T_minus_D_i), F_star+i);
        
        ComplexType Z_i; RealType b_L_D_i;
        sycomore::simd::load_aligned(Z+i, Z_i);
        sycomore::simd::load_aligned(b_L_D+i, b_L_D_i);
        sycomore::simd::store_aligned(Z_i * simd::exp(-b_L_D_i), Z+i);
    }
}

template<INSTRUCTION_SET_TYPE InstructionSet>
void
diffusion_3d_d(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using RealBatch = simd::Batch<Real, InstructionSet>;
    using ComplexBatch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % ComplexBatch::size;
    
    diffusion_3d_w<RealBatch, ComplexBatch>(
        b_L_D, b_T_plus_D, b_T_minus_D, 
        F, F_star, Z, 0, simd_end, ComplexBatch::size);
    diffusion_3d_w<Real, Complex>(
        b_L_D, b_T_plus_D, b_T_minus_D, 
        F, F_star, Z, simd_end, states_count, 1);
}

template<typename ValueType>
void diffusion_3d_b_w(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D, 
    std::size_t begin, std::size_t end, std::size_t step)
{
    for(int i=begin; i<end; i+=step)
    {
        ValueType k_m_i, k_n_i;
        sycomore::simd::load_aligned(k_m+i, k_m_i);
        sycomore::simd::load_aligned(k_n+i, k_n_i);
        
        auto const b_L = k_m_i * k_n_i * tau;
        ValueType b_L_D_i;
        sycomore::simd::load_aligned(b_L_D+i, b_L_D_i);
        sycomore::simd::store_aligned(b_L_D_i+b_L*D_mn, b_L_D+i);
        
        auto const b_T_plus = 
            b_L + delta_k_product_term + 0.5 * tau * (
                k_m_i*delta_k_n + k_n_i*delta_k_m);
        ValueType b_T_plus_D_i;
        sycomore::simd::load_aligned(b_T_plus_D+i, b_T_plus_D_i);
        sycomore::simd::store_aligned(b_T_plus_D_i+b_T_plus*D_mn, b_T_plus_D+i);
        
        auto const b_T_minus = 
            b_L + delta_k_product_term + 0.5 * tau * (
                -k_m_i*delta_k_n + -k_n_i*delta_k_m);
        ValueType b_T_minus_D_i;
        sycomore::simd::load_aligned(b_T_minus_D+i, b_T_minus_D_i);
        sycomore::simd::store_aligned(b_T_minus_D_i+b_T_minus*D_mn, b_T_minus_D+i);
    }
}

template<INSTRUCTION_SET_TYPE InstructionSet>
void
diffusion_3d_b_d(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D,
    unsigned int states_count)
{
    using Batch = simd::Batch<Real, InstructionSet>;
    auto const simd_end = states_count - states_count % Batch::size;
    
    diffusion_3d_b_w<Batch>(
        k_m, k_n, delta_k_m, delta_k_n, delta_k_product_term, tau, D_mn, 
        b_L_D, b_T_plus_D, b_T_minus_D,
        0, simd_end, Batch::size);
    diffusion_3d_b_w<Real>(
        k_m, k_n, delta_k_m, delta_k_n, delta_k_product_term, tau, D_mn, 
        b_L_D, b_T_plus_D, b_T_minus_D,
        simd_end, states_count, 1);
}

/*******************************************************************************
 *                           Off-resonance operator                            *
 ******************************************************************************/

template<typename ValueType>
void off_resonance_w(
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

template<INSTRUCTION_SET_TYPE InstructionSet>
void
off_resonance_d(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using Batch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % Batch::size;
    
    off_resonance_w<Batch>(phi, F, F_star, Z, 0, simd_end, Batch::size);
    off_resonance_w<Complex>(phi, F, F_star, Z, simd_end, states_count, 1);
}

/*******************************************************************************
 *                            Bulk motion operator                             *
 ******************************************************************************/

template<typename RealType, typename ComplexType>
void bulk_motion_w(
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

template<INSTRUCTION_SET_TYPE InstructionSet>
void
bulk_motion_d(
    Real delta_k, Real v, Real tau, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{    
    using RealBatch = simd::Batch<Real, InstructionSet>;
    using ComplexBatch = simd::Batch<Complex, InstructionSet>;
    auto const simd_end = states_count - states_count % ComplexBatch::size;
    
    bulk_motion_w<RealBatch, ComplexBatch>(
        delta_k, v, tau, k, F, F_star, Z, 0, simd_end, ComplexBatch::size);
    bulk_motion_w<Real, Complex>(
        delta_k, v, tau, k, F, F_star, Z, simd_end, states_count, 1);
}

}

}

}

#endif // _c90ddae6_b878_4da1_b3a4_36c025c09ae1
