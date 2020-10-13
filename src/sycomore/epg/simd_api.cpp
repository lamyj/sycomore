#include "simd_api.h"

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

/*******************************************************************************
 *                                Pulse operator                               *
 ******************************************************************************/

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
