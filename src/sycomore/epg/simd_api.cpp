#include "simd_api.h"

#include <utility>
#include <vector>

#include "sycomore/epg/pool_storage.h"
#include "sycomore/simd.h"
#include "sycomore/sycomore.h"

namespace sycomore 
{

namespace epg
{

namespace simd_api
{

/*******************************************************************************
 *                                Pulse operators                              *
 ******************************************************************************/

template<>
void
apply_pulse_single_pool_d<unsupported>(
    std::vector<Complex> const & T,
    pool_storage::SinglePool & storage, unsigned int states_count)
{
    apply_pulse_single_pool_w<Complex>(
        T,
        storage.F.data(), storage.F_star.data(), storage.Z.data(),
        0, states_count, 1);
}

template<>
void
apply_pulse_exchange_d<unsupported>(
    std::vector<Complex> const & T,
    pool_storage::Exchange & storage,
    unsigned int states_count)
{
    apply_pulse_exchange_w<Complex>(
        T,
        storage.F_a.data(), storage.F_star_a.data(), storage.Z_a.data(),
        storage.F_b.data(), storage.F_star_b.data(), storage.Z_b.data(),
        0, states_count, 1);
}

template<>
void
apply_pulse_magnetization_transfer_d<unsupported>(
    std::vector<Complex> const & T,
    pool_storage::MagnetizationTransfer & storage,
    unsigned int states_count)
{
    apply_pulse_magnetization_transfer_w<Complex>(
        T,
        storage.F.data(), storage.F_star.data(), storage.Z_a.data(),
        storage.Z_b.data(), 0, states_count, 1);
}

/*******************************************************************************
 *                             Relaxation operators                            *
 ******************************************************************************/

template<>
void
relaxation_single_pool_d<unsupported>(
    std::pair<Real, Real> const & E,
    pool_storage::SinglePool & storage,
    unsigned int states_count)
{
    // Pass 2*states_count as we are getting reinterpreted real arrays
    relaxation_single_pool_w<Complex>(
        E,
        storage.F.data(), storage.F_star.data(), storage.Z.data(),
        0, 2*states_count, 1);
}

/*******************************************************************************
 *                             Diffusion operator                              *
 ******************************************************************************/

template<>
void
diffusion_d<unsupported>(
    Real delta_k, Real tau, Real D, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    diffusion_w<Real, Complex>(
        delta_k, tau, D, k, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                           3D diffusion operator                             *
 ******************************************************************************/

template<>
void
diffusion_3d_b_d<unsupported>(
    Real const * k_m, Real const * k_n, Real delta_k_m, Real delta_k_n, 
    Real delta_k_product_term, Real tau, Real D_mn,
    Real * b_L_D, Real * b_T_plus_D, Real * b_T_minus_D,
    unsigned int states_count)
{
    diffusion_3d_b_w<Real>(
        k_m, k_n, delta_k_m, delta_k_n, delta_k_product_term, tau, D_mn, 
        b_L_D, b_T_plus_D, b_T_minus_D,
        0, states_count, 1);
}

template<>
void
diffusion_3d_d<unsupported>(
    Real const * b_L_D, Real const * b_T_plus_D, Real const * b_T_minus_D, 
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    diffusion_3d_w<Real, Complex>(
        b_L_D, b_T_plus_D, b_T_minus_D, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                           Off-resonance operator                            *
 ******************************************************************************/

template<>
void
off_resonance_d<unsupported>(
    std::pair<Complex, Complex> const & phi,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    off_resonance_w<Complex>(phi, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                            Bulk motion operator                             *
 ******************************************************************************/

template<>
void
bulk_motion_d<unsupported>(
    Real delta_k, Real v, Real tau, Real const * k,
    Complex * F, Complex * F_star, Complex * Z, unsigned int states_count)
{
    bulk_motion_w<Real, Complex>(
        delta_k, v, tau, k, F, F_star, Z, 0, states_count, 1);
}

/*******************************************************************************
 *                          Function table and set-up                          *
 ******************************************************************************/

decltype(&apply_pulse_single_pool_d<unsupported>)
    apply_pulse_single_pool = nullptr;
decltype(&apply_pulse_exchange_d<unsupported>) apply_pulse_exchange = nullptr;
decltype(&apply_pulse_magnetization_transfer_d<unsupported>)
    apply_pulse_magnetization_transfer = nullptr;
decltype(&relaxation_single_pool_d<unsupported>) relaxation_single_pool = nullptr;
decltype(&diffusion_d<unsupported>) diffusion = nullptr;
decltype(&diffusion_3d_b_d<unsupported>) diffusion_3d_b = nullptr;
decltype(&diffusion_3d_d<unsupported>) diffusion_3d = nullptr;
decltype(&off_resonance_d<unsupported>) off_resonance = nullptr;
decltype(&bulk_motion_d<unsupported>) bulk_motion = nullptr;

void set_api(int instruction_set)
{
    SYCOMORE_SET_API_FUNCTION(apply_pulse_single_pool)
    SYCOMORE_SET_API_FUNCTION(apply_pulse_exchange)
    SYCOMORE_SET_API_FUNCTION(apply_pulse_magnetization_transfer)
    SYCOMORE_SET_API_FUNCTION(relaxation_single_pool)
    SYCOMORE_SET_API_FUNCTION(diffusion)
    SYCOMORE_SET_API_FUNCTION(diffusion_3d_b)
    SYCOMORE_SET_API_FUNCTION(diffusion_3d)
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
