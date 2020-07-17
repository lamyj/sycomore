#include "Regular.h"

#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/Array.h"
#include "sycomore/epg/operators.h"
#include "sycomore/Grid.h"
#include "sycomore/GridScanner.h"
#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

Regular
::Regular(
    Species const & species, Magnetization const & initial_magnetization, 
    unsigned int initial_size, 
    Quantity const & unit_gradient_area, double gradient_tolerance)
: species(species),
    _F(initial_size, 0), _F_star(initial_size, 0), _Z(initial_size, 0),
    _unit_gradient_area(unit_gradient_area), _gradient_tolerance(gradient_tolerance)
{
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_F[0] = std::sqrt(2)*magnetization.p;
    this->_F_star[0] = std::sqrt(2)*magnetization.m;
    this->_Z[0] = magnetization.z;
    
    this->_states_count = 1;
}

std::size_t const 
Regular
::states_count() const
{
    return this->_states_count;
}

std::vector<Complex>
Regular
::state(std::size_t order) const
{
    return {this->_F[order], this->_F_star[order], this->_Z[order]};
}

std::vector<Complex>
Regular
::states() const
{
    std::vector<Complex> result(3*this->_states_count);
    for(unsigned int order=0; order<this->_states_count; ++order)
    {
        result[3*order+0] = this->_F[order];
        result[3*order+1] = this->_F_star[order];
        result[3*order+2] = this->_Z[order];
    }
    return result;
}

Complex const &
Regular
::echo() const
{
    return this->_F[0];
}

using simd_type = xsimd::simd_type<Complex>;
constexpr std::size_t const simd_width = simd_type::size;

void
Regular
::apply_pulse(Quantity angle, Quantity phase)
{
    auto const T = operators::pulse(angle.magnitude, phase.magnitude);
    
    std::size_t const simd_size = this->_states_count - this->_states_count % simd_width;
    
    for(std::size_t order = 0; order < simd_size; order += simd_width)
    {
        auto const F = xsimd::load_aligned(&this->_F[order]);
        auto const F_star = xsimd::load_aligned(&this->_F_star[order]);
        auto const Z = xsimd::load_aligned(&this->_Z[order]);
    
        auto const F_new =      T[3*0+0] * F + T[3*0+1] * F_star + T[3*0+2] * Z;
        auto const F_star_new = T[3*1+0] * F + T[3*1+1] * F_star + T[3*1+2] * Z;
        auto const Z_new =      T[3*2+0] * F + T[3*2+1] * F_star + T[3*2+2] * Z;
    
        xsimd::store_aligned(&this->_F[order], F_new);
        xsimd::store_aligned(&this->_F_star[order], F_star_new);
        xsimd::store_aligned(&this->_Z[order], Z_new);
    }
    
    for(std::size_t order = simd_size; order < this->_states_count; ++order)
    {
        auto const & F = _F[order];
        auto const & F_star = _F_star[order];
        auto const & Z = _Z[order];
    
        auto const F_new =      T[3*0+0] * F + T[3*0+1] * F_star + T[3*0+2] * Z;
        auto const F_star_new = T[3*1+0] * F + T[3*1+1] * F_star + T[3*1+2] * Z;
        auto const Z_new =      T[3*2+0] * F + T[3*2+1] * F_star + T[3*2+2] * Z;
    
        this->_F[order] = F_new;
        this->_F_star[order] = F_star_new;
        this->_Z[order] = Z_new;
    }
}

void
Regular
::apply_time_interval(Quantity const & duration, Quantity const & gradient)
{
    // Note that since E does not depend on k, the E and S operators commute
    // and that E and D(k) also commute as they are diagonal matrices. The
    // only effect will be the relative order of D and S.
    // Since the diffusion operator relies on the "start" state k_1, we need
    // to apply the gradient operator after the diffusion operator. Otherwise
    // states would be dephased by D(k+Δk, Δk) instead of D(k, Δk)
    
    this->relaxation(duration);
    this->diffusion(duration, gradient);
    this->bulk_motion(duration, gradient);
    if(duration.magnitude != 0 && gradient.magnitude != 0)
    {
        if(this->_unit_gradient_area.magnitude != 0)
        {
            this->shift(duration, gradient);
        }
        else
        {
            this->shift();
        }
    }
    if(
        duration.magnitude != 0 
        && (this->delta_omega.magnitude != 0 || species.get_delta_omega().magnitude != 0))
    {
        this->off_resonance(duration);
    }
}

void
Regular
::apply_time_interval(TimeInterval const & interval)
{
    this->apply_time_interval(
        interval.get_duration(), interval.get_gradient_amplitude()[0]);
}

void
Regular
::shift()
{
    this->_shift(1);
}

void
Regular
::shift(Quantity const & duration, Quantity const & gradient)
{
    auto const area = duration*gradient;
    auto const epsilon = this->_gradient_tolerance*this->_unit_gradient_area.magnitude;
    auto const remainder = std::remainder(
        area.magnitude, this->_unit_gradient_area.magnitude);
    
    if(std::abs(remainder) >= epsilon)
    {
        throw std::runtime_error(
            "Gradient is not a interger multiple of unit gradient");
    }
    
    int n = std::round(area/this->_unit_gradient_area);
    
    this->_shift(n);
}

void
Regular
::relaxation(Quantity const & duration)
{
    if(this->species.get_R1().magnitude == 0 && this->species.get_R2().magnitude == 0)
    {
        return;
    }
    
    auto const E = operators::relaxation(
        this->species.get_R1().magnitude, this->species.get_R2().magnitude, 
        duration.magnitude);
    // SIMD operations on Real data are faster than those on Complex data.
    using simd_type = xsimd::simd_type<Real>;
    constexpr std::size_t const simd_width = simd_type::size;
    auto const size = 2*_states_count;
    std::size_t const simd_size = size - size % simd_width;
    
    auto F = reinterpret_cast<Real*>(this->_F.data());
    auto F_star = reinterpret_cast<Real*>(this->_F_star.data());
    auto Z = reinterpret_cast<Real*>(this->_Z.data());
    
    for(std::size_t i = 0; i < simd_size; i += simd_width)
    {
        xsimd::store_aligned(&F[i], xsimd::load_aligned(&F[i])*E.second);
        xsimd::store_aligned(&F_star[i], xsimd::load_aligned(&F_star[i])*E.second);
        xsimd::store_aligned(&Z[i], xsimd::load_aligned(&Z[i])*E.first);
    }
    
    for(std::size_t i = simd_size; i < size; ++i)
    {
        F[i] *= E.second;
        F_star[i] *= E.second;
        Z[i] *= E.first;
    }
    
    this->_Z[0] += 1.-E.first; // WARNING: assumes M0=1
}

void
Regular
::diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(this->species.get_D()[0].magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = (sycomore::gamma*gradient*duration).magnitude;
    if(delta_k == 0)
    {
        return;
    }
    
    // NOTE: the diffusion operator is real, but running it on the real 
    // components of the complex data as for the relaxation operator requires
    // computing the operator for each order twice and the tweak becomes slower
    // than working directly on the complex data.
    
    std::vector<Real> k_array(this->_states_count);
    for(std::size_t i=0; i<k_array.size(); ++i)
    {
        k_array[i] = delta_k*i;
    }
    
    auto const & tau = duration.magnitude;
    auto const & D = species.get_D()[0].magnitude;
    
    std::size_t const simd_size = this->_states_count - this->_states_count % simd_width;
    // #pragma omp parallel for
    for(std::size_t i = 0; i < simd_size; i += simd_width)
    {
        auto const k = xsimd::load_aligned(&k_array[i]);
    
        auto const F = xsimd::load_aligned(&this->_F[i]);
        auto const b_T_plus = tau*(xsimd::pow(k+delta_k/2, 2) + std::pow(delta_k, 2) / 12);
        auto const D_T_plus = xsimd::exp(-b_T_plus*D);
        xsimd::store_aligned(&this->_F[i], F*D_T_plus);
    
        auto const F_star = xsimd::load_aligned(&this->_F_star[i]);
        auto const b_T_minus = tau*(xsimd::pow(-k+delta_k/2, 2) + std::pow(delta_k, 2) / 12);
        auto const D_T_minus = xsimd::exp(-b_T_minus*D);
        xsimd::store_aligned(&this->_F_star[i], F_star*D_T_minus);
    
        auto const Z = xsimd::load_aligned(&this->_Z[i]);
        auto const b_L = xsimd::pow(k, 2) * tau;
        auto const D_L = xsimd::exp(-b_L*D);
        xsimd::store_aligned(&this->_Z[i], Z*D_L);
    }
    
    for(std::size_t i = simd_size; i < _states_count; ++i)
    {
        auto const & k = k_array[i];
    
        auto const b_T_plus = tau*(std::pow(k+delta_k/2, 2.) + std::pow(delta_k, 2.) / 12);
        auto const D_T_plus = std::exp(-b_T_plus*D);
        this->_F[i] *= D_T_plus;
    
        auto const b_T_minus = tau*(std::pow(-k+delta_k/2, 2.) + std::pow(delta_k, 2.) / 12);
        auto const D_T_minus = std::exp(-b_T_minus*D);
        this->_F_star[i] *= D_T_minus;
    
        auto const b_L = std::pow(k, 2) * tau;
        auto const D_L = std::exp(-b_L*D);
        this->_Z[i] *= D_L;
    }
}

void
Regular
::off_resonance(Quantity const & duration)
{
    auto const angle = 
        duration * 2*M_PI*units::rad 
        * (this->delta_omega+this->species.get_delta_omega());
    if(angle.magnitude != 0)
    {
        auto const rotations = operators::phase_accumulation(angle.magnitude);
        
        for(int order=0; order<this->_states_count; ++order)
        {
            this->_F[order] *= rotations.first;
            this->_F_star[order] *= rotations.second;
            // Z̃ states are unaffected
        }
    }
}

void
Regular
::bulk_motion(Quantity const & duration, Quantity const & gradient)
{
    if(velocity.magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = (sycomore::gamma*gradient*duration).magnitude;
    if(delta_k == 0)
    {
        return;
    }
    
    #pragma omp parallel for schedule(static)
    for(int order=0; order<this->_states_count; ++order)
    {
        auto const k = order*delta_k;
        auto const J = operators::bulk_motion(
            this->velocity.magnitude, duration.magnitude, k, delta_k);
        this->_F[order] *= std::get<0>(J);
        this->_F_star[order] *= std::get<1>(J);
        this->_Z[order] *= std::get<2>(J);
    }
}

Quantity const &
Regular
::unit_gradient_area() const
{
    return this->_unit_gradient_area;
}

double
Regular
::gradient_tolerance() const
{
    return this->_gradient_tolerance;
}

void
Regular
::_shift(int n)
{
    if(std::abs(n)>1)
    {
        for(int i=0; i<std::abs(n); ++i)
        {
            this->_shift(n>0?+1:-1);
        }
    }
    else if(n==1 || n==-1)
    {
        if(this->_states_count >= this->_F.size())
        {
            this->_F.resize(this->_F.size()+100, 0);
            this->_F_star.resize(this->_F_star.size()+100, 0);
            this->_Z.resize(this->_Z.size()+100, 0);
        }
        
        if(n == +1)
        {
            // Shift positive F̃ states right
            std::copy_backward(
                this->_F.begin(), this->_F.begin()+this->_states_count, 
                this->_F.begin()+this->_states_count+1);
            
            // Shift negative F̃^* states left
            std::copy(
                this->_F_star.begin()+1, 
                this->_F_star.begin()+this->_states_count+1, 
                this->_F_star.begin());
            
            // Update F̃_{+0} using F̃^*_{-0}
            this->_F[0] = std::conj(this->_F_star[0]);
        }
        else
        {
            // Shift negative F̃^* states right
            std::copy_backward(
                this->_F_star.begin(), 
                this->_F_star.begin()+this->_states_count, 
                this->_F_star.begin()+this->_states_count+1);
            
            // Shift positive F̃ states left
            std::copy(
                this->_F.begin()+1, 
                this->_F.begin()+this->_states_count+1, 
                this->_F.begin());
            
            // Update F̃^*_{-0} using F̃_{+0}
            this->_F_star[0] = std::conj(this->_F[0]);
        }
        
        ++this->_states_count;
        
        // Remove empty states with high order.
        while(
            this->_F[this->_states_count-1] == 0. 
            && this->_F_star[this->_states_count-1] == 0.
            && this->_Z[this->_states_count-1] == 0.)
        {
            --this->_states_count;
        }
    }
    else
    { 
        // n==0, nothing to do.
    }
}

}

}
