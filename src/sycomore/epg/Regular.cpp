#include "Regular.h"

#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/epg/operators.h"
#include "sycomore/epg/simd_api.h"
#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/TimeInterval.h"
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
    _unit_gradient_area(unit_gradient_area), 
    _gradient_tolerance(gradient_tolerance)
{
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_F[0] = std::sqrt(2)*magnetization.p;
    this->_F_star[0] = std::sqrt(2)*magnetization.m;
    this->_Z[0] = magnetization.z;
    this->_M_z_eq = magnetization.z;
    
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

void
Regular
::apply_pulse(Quantity angle, Quantity phase)
{
    simd_api::apply_pulse(
        operators::pulse(angle.magnitude, phase.magnitude), 
        this->_F.data(), this->_F_star.data(), this->_Z.data(), 
        this->_states_count);
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
        && (
            this->delta_omega.magnitude != 0 
            || species.get_delta_omega().magnitude != 0))
    {
        this->off_resonance(duration);
    }
    
    if(threshold > 0)
    {
        // Remove low-populated states with high order.
        auto const threshold_squared = std::pow(this->threshold, 2);
        bool done = false;
        while(this->_states_count>1 && !done)
        {
            auto const magnitude_squared =
                std::pow(std::abs(this->_F[this->_states_count-1]), 2)
                +std::pow(std::abs(this->_F_star[this->_states_count-1]), 2)
                +std::pow(std::abs(this->_Z[this->_states_count-1]), 2);
            
            if(magnitude_squared > threshold_squared)
            {
                done = true;
            }
            else
            {
                --this->_states_count;
            }
        }
    }
    else
    {
        // Remove empty states with high order.
        while(
            this->_states_count > 1
            && this->_F[this->_states_count-1] == 0. 
            && this->_F_star[this->_states_count-1] == 0.
            && this->_Z[this->_states_count-1] == 0.)
        {
            --this->_states_count;
        }
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
    auto const epsilon = 
        this->_gradient_tolerance*this->_unit_gradient_area.magnitude;
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
    if(
        this->species.get_R1().magnitude == 0 
        && this->species.get_R2().magnitude == 0)
    {
        return;
    }
    
    auto const E = operators::relaxation(
        this->species.get_R1().magnitude, this->species.get_R2().magnitude, 
        duration.magnitude);
    
    simd_api::relaxation(
        E,
        reinterpret_cast<Real*>(this->_F.data()),
        reinterpret_cast<Real*>(this->_F_star.data()),
        reinterpret_cast<Real*>(this->_Z.data()),
        this->_states_count);
    
    this->_Z[0] += this->_M_z_eq*(1.-E.first);
}

void
Regular
::diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(this->species.get_D()[0].magnitude == 0)
    {
        return;
    }
    
    auto const area = duration.magnitude*gradient.magnitude;
    if(area == 0)
    {
        return;
    }
    
    auto const unit_gradient_area = this->_unit_gradient_area.magnitude;
    if(unit_gradient_area == 0)
    {
        throw std::runtime_error(
            "Cannot compute diffusion without unit gradient area");
    }
    
    auto const remainder = std::remainder(area, unit_gradient_area);
    if(std::abs(remainder) >= this->_gradient_tolerance*unit_gradient_area)
    {
        throw std::runtime_error(
            "Gradient is not a interger multiple of unit gradient");
    }
    
    auto const delta_k = sycomore::gamma.magnitude*area;
    
    this->_cache.update_diffusion(this->_states_count, unit_gradient_area);
    
    auto const & tau = duration.magnitude;
    auto const & D = species.get_D()[0].magnitude;
    
    simd_api::diffusion(
        delta_k, tau, D, this->_cache.k.data(),
        this->_F.data(), this->_F_star.data(), this->_Z.data(),
        this->_states_count);
}

void
Regular
::off_resonance(Quantity const & duration)
{
    auto const angle = 
        duration.magnitude * 2*M_PI*units::rad 
        * (this->delta_omega.magnitude+this->species.get_delta_omega().magnitude);
    if(angle != 0)
    {
        auto const rotations = operators::phase_accumulation(angle);
        simd_api::off_resonance(
            rotations,
            this->_F.data(), this->_F_star.data(), this->_Z.data(),
            this->_states_count);
    }
}

void
Regular
::bulk_motion(Quantity const & duration, Quantity const & gradient)
{
    if(this->velocity.magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = (
        sycomore::gamma.magnitude*gradient.magnitude*duration.magnitude);
    if(delta_k == 0)
    {
        return;
    }
    
    std::vector<Real, xsimd::aligned_allocator<Real, 64>> k(
        this->_states_count);
    for(std::size_t i=0; i<k.size(); ++i)
    {
        k[i] = delta_k*i;
    }
    
    simd_api::bulk_motion(
        delta_k, this->velocity.magnitude, duration.magnitude, k.data(),
        this->_F.data(), this->_F_star.data(), this->_Z.data(),
        this->_states_count);
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
            // Shift positive F states right
            std::copy_backward(
                this->_F.begin(), this->_F.begin()+this->_states_count, 
                this->_F.begin()+this->_states_count+1);
            
            // Shift negative F* states left
            std::copy(
                this->_F_star.begin()+1, 
                this->_F_star.begin()+this->_states_count+1, 
                this->_F_star.begin());
            
            // Update extremal states: F_{+0} using F*_{-0}, F*_{-max+1}=0
            this->_F[0] = std::conj(this->_F_star[0]);
            this->_F_star[this->_states_count] = 0;
        }
        else
        {
            // Shift negative F* states right
            std::copy_backward(
                this->_F_star.begin(), 
                this->_F_star.begin()+this->_states_count, 
                this->_F_star.begin()+this->_states_count+1);
            
            // Shift positive F states left
            std::copy(
                this->_F.begin()+1, 
                this->_F.begin()+this->_states_count+1, 
                this->_F.begin());
            
            // Update extremal states: F*_{-0} using F_{+0}, F_{max+1}=0
            this->_F_star[0] = std::conj(this->_F[0]);
            this->_F[this->_states_count] = 0;
        }
        
        ++this->_states_count;
    }
    else
    { 
        // n==0, nothing to do.
    }
}

void
Regular::Cache
::update_diffusion(std::size_t size, Real unit_gradient_area)
{
    this->k.resize(size);
    for(std::size_t order=0; order != size; ++order)
    {
        this->k[order] = order*sycomore::gamma.magnitude*unit_gradient_area;
    }
}

}

}
