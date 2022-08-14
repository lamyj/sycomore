#include "Regular.h"

#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/epg/Base.h"
#include "sycomore/epg/operators.h"
#include "sycomore/epg/pool_model.h"
#include "sycomore/epg/pool_storage.h"
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
: Base(species, initial_magnetization, initial_size),
    _unit_gradient_area(unit_gradient_area), 
    _gradient_tolerance(gradient_tolerance)
{
    this->_states_count = 1;
}

std::size_t 
Regular
::size() const
{
    return this->_states_count;
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
            || this->_model->species.get_delta_omega().magnitude != 0))
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
                std::pow(std::abs(this->_storage->F[this->_states_count-1]), 2)
                +std::pow(std::abs(this->_storage->F_star[this->_states_count-1]), 2)
                +std::pow(std::abs(this->_storage->Z[this->_states_count-1]), 2);
            
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
            && this->_storage->F[this->_states_count-1] == 0. 
            && this->_storage->F_star[this->_states_count-1] == 0.
            && this->_storage->Z[this->_states_count-1] == 0.)
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
::diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(this->_model->species.get_D()[0].magnitude == 0)
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
    auto const & D = this->_model->species.get_D()[0].magnitude;
    
    simd_api::diffusion(
        delta_k, tau, D, this->_cache.k.data(),
        this->_storage->F.data(),
        this->_storage->F_star.data(),
        this->_storage->Z.data(),
        this->_states_count);
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
        this->_storage->F.data(),
        this->_storage->F_star.data(),
        this->_storage->Z.data(),
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
        auto & [F, F_star, Z] = *this->_storage;
        if(this->_states_count >= F.size())
        {
            F.resize(F.size()+100, 0);
            F_star.resize(F_star.size()+100, 0);
            Z.resize(Z.size()+100, 0);
        }
        
        if(n == +1)
        {
            // Shift positive F states right
            std::copy_backward(
                F.begin(), F.begin()+this->_states_count, 
                F.begin()+this->_states_count+1);
            
            // Shift negative F* states left
            std::copy(
                F_star.begin()+1, F_star.begin()+this->_states_count+1, 
                F_star.begin());
            
            // Update extremal states: F_{+0} using F*_{-0}, F*_{-max+1}=0
            F[0] = std::conj(F_star[0]);
            F_star[this->_states_count] = 0;
        }
        else
        {
            // Shift negative F* states right
            std::copy_backward(
                F_star.begin(), F_star.begin()+this->_states_count, 
                F_star.begin()+this->_states_count+1);
            
            // Shift positive F states left
            std::copy(F.begin()+1, F.begin()+this->_states_count+1, F.begin());
            
            // Update extremal states: F*_{-0} using F_{+0}, F_{max+1}=0
            F_star[0] = std::conj(F[0]);
            F[this->_states_count] = 0;
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
