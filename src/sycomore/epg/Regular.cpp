#include "Regular.h"

#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/Array.h"
#include "sycomore/epg/Base.h"
#include "sycomore/epg/operators.h"
#include "sycomore/epg/simd_api.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

Regular
::Regular(
    Species const & species, Vector3R const & initial_magnetization, 
    unsigned int initial_size, 
    Quantity const & unit_dephasing, double gradient_tolerance)
: Base(species, initial_magnetization, initial_size),
    _states_count(1), _unit_dephasing(unit_dephasing), 
    _gradient_tolerance(gradient_tolerance)
{
    if(this->_unit_dephasing.dimensions != GradientDephasing)
    {
        this->_unit_dephasing *= sycomore::gamma;
    }
}

Regular
::Regular(
    Species const & species_a, Species const & species_b,
    Vector3R const & M0_a, Vector3R const & M0_b,
    Quantity const & k_a, Quantity const & delta_b,
    unsigned int initial_size, Quantity const & unit_dephasing,
    double gradient_tolerance)
: Base(species_a, species_b, M0_a, M0_b, k_a, delta_b, initial_size),
    _states_count(1), _unit_dephasing(unit_dephasing), 
    _gradient_tolerance(gradient_tolerance)
{
    if(this->_unit_dephasing.dimensions != GradientDephasing)
    {
        this->_unit_dephasing *= sycomore::gamma;
    }
}

Regular
::Regular(
    Species const & species_a, Quantity const & R1_b_or_T1_b,
    Vector3R const & M0_a, Vector3R const & M0_b,
    Quantity const & k_a,
    unsigned int initial_size, Quantity const & unit_dephasing,
    double gradient_tolerance)
: Base(species_a, R1_b_or_T1_b, M0_a, M0_b, k_a, initial_size),
    _states_count(1), _unit_dephasing(unit_dephasing), 
    _gradient_tolerance(gradient_tolerance)
{
    if(this->_unit_dephasing.dimensions != GradientDephasing)
    {
        this->_unit_dephasing *= sycomore::gamma;
    }
}

std::size_t 
Regular
::size() const
{
    return this->_states_count;
}

TensorQ<1>
Regular
::orders() const
{
    TensorQ<1> result(TensorQ<1>::shape_type{this->size()});
    auto const factor = 
        (this->_unit_dephasing.magnitude != 0)
        ? this->_unit_dephasing : Quantity(1.);
    for(std::size_t i=0; i<result.size(); ++i)
    {
        result[i] = i * factor;
    }
    return result;
}

ArrayC
Regular
::state(Order const & order) const
{
    std::size_t position;
    if(order.dimensions == Dimensions())
    {
        position = std::size_t(order.magnitude);
    }
    else
    {
        auto const epsilon = 
            this->_gradient_tolerance*this->_unit_dephasing.magnitude;
        auto const remainder = std::remainder(
            order.magnitude, this->_unit_dephasing.magnitude);
        if(std::abs(remainder) >= epsilon)
        {
            throw std::runtime_error(
                "Dephasing is not a integer multiple of unit dephasing");
        }
        position = std::round(order/this->_unit_dephasing);
    }
    
    return this->state(position);
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
        if(this->_unit_dephasing.magnitude != 0)
        {
            this->shift(duration, gradient);
        }
        else
        {
            this->shift();
        }
    }
    this->off_resonance(duration);
    
    this->_elapsed += duration.magnitude;
    
    // Remove low-populated states with high order.
    auto const threshold_squared = std::pow(this->threshold, 2);
    
    bool done = false;
    while(this->_states_count > 1 && !done)
    {
        Real max_magnitude_squared = 0.;
        for(std::size_t pool=0; pool<this->_model.pools; ++pool)
        {
            using std::pow; using std::abs;
            auto const magnitude_squared = 
                pow(abs(this->_model.F[pool][this->_states_count-1]), 2)
                +pow(abs(this->_model.F_star[pool][this->_states_count-1]), 2)
                +pow(abs(this->_model.Z[pool][this->_states_count-1]), 2);
            max_magnitude_squared = std::max(
                max_magnitude_squared, magnitude_squared);
        }
        
        if(max_magnitude_squared > threshold_squared)
        {
            done = true;
        }
        else
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
        interval.duration(), interval.gradient_amplitude()[0]);
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
    auto const dephasing = sycomore::gamma*duration*gradient;
    auto const epsilon = 
        this->_gradient_tolerance*this->_unit_dephasing.magnitude;
    auto const remainder = std::remainder(
        dephasing.magnitude, this->_unit_dephasing.magnitude);
    
    if(std::abs(remainder) >= epsilon)
    {
        throw std::runtime_error(
            "Dephasing is not a integer multiple of unit dephasing");
    }
    
    int n = std::round(dephasing/this->_unit_dephasing);
    
    this->_shift(n);
}

void
Regular
::diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(
        std::all_of(
            this->_model.species.begin(), this->_model.species.end(),
            [](Species const & s) { return s.D().unchecked(0, 0).magnitude == 0; }))
    {
        return;
    }
    
    auto const dephasing =
        sycomore::gamma.magnitude * duration.magnitude * gradient.magnitude;
    if(dephasing == 0)
    {
        return;
    }
    
    auto const unit_dephasing = this->_unit_dephasing.magnitude;
    if(unit_dephasing == 0)
    {
        throw std::runtime_error(
            "Cannot compute diffusion without unit dephasing");
    }
    
    auto const remainder = std::remainder(dephasing, unit_dephasing);
    if(std::abs(remainder) >= this->_gradient_tolerance*unit_dephasing)
    {
        throw std::runtime_error(
            "Gradient is not a interger multiple of unit gradient");
    }
    
    auto const delta_k = dephasing;
    
    this->_cache.update_diffusion(this->size(), unit_dephasing);
    
    auto const & tau = duration.magnitude;
    
    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
    {
        auto const & D = this->_model.species[pool].D().unchecked(0, 0).magnitude;
        if(D == 0.)
        {
            continue;
        }
        
        simd_api::diffusion(
            delta_k, tau, D, this->_cache.k.data(),
            this->_model.F[pool], this->_model.F_star[pool], this->_model.Z[pool],
            this->size());
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
    
    std::vector<Real, xsimd::aligned_allocator<Real, 64>> k(this->size());
    for(std::size_t i=0; i<k.size(); ++i)
    {
        k[i] = delta_k*i;
    }
    
    simd_api::bulk_motion(
        delta_k, this->velocity.magnitude, duration.magnitude, k.data(),
        this->_model, this->size());
}

Quantity const &
Regular
::unit_dephasing() const
{
    return this->_unit_dephasing;
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
        auto const size = this->size();
        for(std::size_t pool=0; pool<this->_model.pools; ++pool)
        {
            auto & F = this->_model.F[pool];
            auto & F_star = this->_model.F_star[pool];
            auto & Z = this->_model.Z[pool];
            if(size >= F.size())
            {
                F.resize(F.size()+100, 0);
                F_star.resize(F_star.size()+100, 0);
                Z.resize(Z.size()+100, 0);
            }
            
            if(n == +1)
            {
                // Shift positive F states right
                std::copy_backward(F.begin(), F.begin()+size, F.begin()+size+1);
                
                // Shift negative F* states left
                std::copy(
                    F_star.begin()+1, F_star.begin()+size+1, F_star.begin());
                
                // Update extremal states: F_{+0} using F*_{-0}, F*_{-max+1}=0
                F[0] = std::conj(F_star[0]);
                F_star[size] = 0;
            }
            else
            {
                // Shift negative F* states right
                std::copy_backward(
                    F_star.begin(), F_star.begin()+size, F_star.begin()+size+1);
                
                // Shift positive F states left
                std::copy(F.begin()+1, F.begin()+size+1, F.begin());
                
                // Update extremal states: F*_{-0} using F_{+0}, F_{max+1}=0
                F_star[0] = std::conj(F[0]);
                F[size] = 0;
            }
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
::update_diffusion(std::size_t size, Real unit_dephasing)
{
    this->k.resize(size);
    for(std::size_t order=0; order != size; ++order)
    {
        this->k[order] = order*unit_dephasing;
    }
}

}

}
