#include "Discrete.h"

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "sycomore/epg/operators.h"
#include "sycomore/epg/robin_hood.h"
#include "sycomore/epg/simd_api.h"
#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

Discrete
::Discrete(
    Species const & species, Magnetization const & initial_magnetization, 
    Quantity bin_width, Real threshold)
: species(species), _bin_width(bin_width), _F(1, 0), _F_star(1, 0), _Z(1, 0), 
    threshold(threshold)
{
    // Store magnetization as lines of F, F*_, Z
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_F[0] = std::sqrt(2)*magnetization.p;
    this->_F_star[0] = std::sqrt(2)*magnetization.m;
    this->_Z[0] = magnetization.z;
    this->_M_z_eq = magnetization.z;
    
    this->_orders.push_back(0);
}

std::size_t
Discrete
::size() const
{
    return this->_orders.size();
}

std::vector<Quantity>
Discrete
::orders() const
{
    std::vector<Quantity> orders(this->_orders.size());
    std::transform(
        this->_orders.begin(), this->_orders.end(), orders.begin(),
        [&](decltype(this->_orders[0]) k){ return k*this->_bin_width; });
    return orders;
}

std::vector<Complex>
Discrete
::state(std::size_t order) const
{
    return {this->_F[order], this->_F_star[order], this->_Z[order]};
}

std::vector<Complex>
Discrete
::state(Quantity const & order) const
{
    using namespace sycomore::units;

    std::size_t const k = std::round(double(order/this->_bin_width));

    auto const it = std::find(this->_orders.begin(), this->_orders.end(), k);
    if(it == this->_orders.end())
    {
        std::ostringstream message;
        message << "No such order: " << order;
        throw std::runtime_error(message.str());
    }

    auto const index = it-this->_orders.begin();
    return this->state(index);
}

std::vector<Complex>
Discrete
::states() const
{
    std::vector<Complex> result(3*this->_orders.size());
    for(unsigned int order=0; order<this->_orders.size(); ++order)
    {
        result[3*order+0] = this->_F[order];
        result[3*order+1] = this->_F_star[order];
        result[3*order+2] = this->_Z[order];
    }
    return result;
}

Complex const &
Discrete
::echo() const
{
    return this->_F[0];
}

void
Discrete
::apply_pulse(Quantity angle, Quantity phase)
{
    simd_api::apply_pulse(
        operators::pulse(angle.magnitude, phase.magnitude), 
        this->_F.data(), this->_F_star.data(), this->_Z.data(), 
        this->_orders.size());
}

void
Discrete
::apply_time_interval(
    Quantity const & duration, Quantity const & gradient, Real threshold)
{
    static bool warning_displayed = false;
    if(!warning_displayed && threshold > 0)
    {
        std::cout 
            << "WARNING: threshold argument is deprecated "
            << "and will be removed from Discrete::apply_time_interval\n";
        warning_displayed = true;
    }
    
    if(duration.magnitude == 0)
    {
        return;
    }
    
    this->relaxation(duration);
    this->diffusion(duration, gradient);
    this->bulk_motion(duration, gradient);
    this->shift(duration, gradient);
    this->off_resonance(duration);
    
    if(threshold == 0)
    {
        threshold = this->threshold;
    }
    if(threshold > 0)
    {
        auto const threshold_squared = std::pow(threshold, 2);
        
        std::size_t destination=1;
        for(std::size_t source=1, end=this->size(); source != end; ++source)
        {
            auto const magnitude_squared =
                std::pow(std::abs(this->_F[source]), 2)
                +std::pow(std::abs(this->_F_star[source]), 2)
                +std::pow(std::abs(this->_Z[source]), 2);
            // Always include the zero order (implicit since we start at 1),
            // include other order if population is above threshold.
            if(magnitude_squared >= threshold_squared)
            {
                if(source != destination)
                {
                    this->_orders[destination] = this->_orders[source];
                    this->_F[destination] = this->_F[source];
                    this->_F_star[destination] = this->_F_star[source];
                    this->_Z[destination] = this->_Z[source];
                }
                ++destination;
            }
        }

        this->_orders.resize(destination);
        this->_F.resize(destination);
        this->_F_star.resize(destination);
        this->_Z.resize(destination);
        
        // No need to update the iterator pointing to the echo magnetization.
    }
}

void
Discrete
::apply_time_interval(TimeInterval const & interval)
{
    this->apply_time_interval(
        interval.get_duration(), interval.get_gradient_amplitude()[0], 0);
}

void
Discrete
::shift(Quantity const & duration, Quantity const & gradient)
{
    // This assumes a constant gradient in the integral: 
    // k(t) = γ ∫_0^t G(t') dt' = γ⋅t⋅G
    long long const delta_k = std::round(
        sycomore::gamma.magnitude*gradient.magnitude*duration.magnitude 
        / this->_bin_width.magnitude);
    
    if(delta_k == 0)
    {
        return;
    }
    
    this->_cache.update_shift(this->size());
    
    for(std::size_t i=0, end=this->size(); i != end; ++i)
    {
        auto const k = this->_orders[i];
        
        auto state = this->_Z[i];
        if(state != 0.)
        {
            this->_cache.Z[this->_cache.get_location(k)] = state;
        }
        
        state = this->_F[i];
        if(state != 0.)
        {
            auto k_F = k+delta_k;
            
            // Depending on whether the new order changed half space, conjugate
            // the state and store it in F* instead of F.
            auto destination = &this->_cache.F;
            if(k_F < 0)
            {
                k_F *= -1;
                destination = &this->_cache.F_star;
                state = std::conj(state);
            }
            (*destination)[this->_cache.get_location(k_F)] = state;
        }
        
        // WARNING: F* state at echo is a duplicate of F state.
        state = this->_F_star[i];
        if(i != 0 && state != 0.)
        {
            // The F* order corresponding to F order k+Δk is -(-k+Δk), i.e. k-Δk
            auto k_F_star = k-delta_k;
            
            // Same as above.
            auto destination = &this->_cache.F_star;
            if(k_F_star < 0)
            {
                k_F_star *= -1;
                destination = &this->_cache.F;
                state = std::conj(state);
            }
            (*destination)[this->_cache.get_location(k_F_star)] = state;
        }
    }
    
    // Update the current orders and states with the new ones.
    this->_cache.orders.resize(this->_cache.locations.size());
    this->_cache.F.resize(this->_cache.locations.size());
    this->_cache.F_star.resize(this->_cache.locations.size());
    this->_cache.Z.resize(this->_cache.locations.size());
    
    // Use swap and not move since we keep the temporary variables between runs.
    std::swap(this->_cache.orders, this->_orders);
    std::swap(this->_cache.F, this->_F);
    std::swap(this->_cache.F_star, this->_F_star);
    std::swap(this->_cache.Z, this->_Z);
    
    // Update the conjugate states of the echo magnetization.
    if(this->_F[0] != 0.)
    {
        this->_F_star[0] = std::conj(this->_F[0]);
    }
    else if(this->_F_star[0] != 0.)
    {
        this->_F[0] = std::conj(this->_F_star[0]);
    }
}

void
Discrete
::relaxation(Quantity const & duration)
{
    if(this->species.get_R1().magnitude == 0 && this->species.get_R2().magnitude == 0)
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
        this->_orders.size());
    
    this->_Z[0] += this->_M_z_eq*(1.-E.first);
}

void
Discrete
::diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(this->species.get_D()[0].magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = (
        sycomore::gamma.magnitude*gradient.magnitude*duration.magnitude);
    if(delta_k == 0)
    {
        return;
    }
    
    this->_cache.update_diffusion(
        this->size(), this->_orders, this->_bin_width.magnitude);
    
    auto const & tau = duration.magnitude;
    auto const & D = species.get_D()[0].magnitude;
    
    simd_api::diffusion(
        delta_k, tau, D, this->_cache.k.data(),
        this->_F.data(), this->_F_star.data(), this->_Z.data(),
        this->_orders.size());
}

void
Discrete
::off_resonance(Quantity const & duration)
{
    auto const angle = 
        duration.magnitude * 2*M_PI
        * (this->delta_omega.magnitude+this->species.get_delta_omega().magnitude);
    if(angle != 0)
    {
        auto const rotations = operators::phase_accumulation(angle);
        simd_api::off_resonance(
            rotations,
            this->_F.data(), this->_F_star.data(), this->_Z.data(),
            this->_orders.size());
    }
}

void
Discrete
::bulk_motion(Quantity const & duration, Quantity const & gradient)
{
    if(this->velocity.magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = (sycomore::gamma*gradient*duration).magnitude;
    if(delta_k == 0)
    {
        return;
    }
    
    std::vector<Real, xsimd::aligned_allocator<Real, 64>> k(this->size());
    for(std::size_t i=0; i<k.size(); ++i)
    {
        k[i] = this->_orders[i] * this->_bin_width.magnitude;
    }
    
    simd_api::bulk_motion(
        delta_k, this->velocity.magnitude, duration.magnitude, k.data(),
        this->_F.data(), this->_F_star.data(), this->_Z.data(),
        this->size());
}

Quantity const & 
Discrete
::bin_width() const
{
    return this->_bin_width;
}

void
Discrete::Cache
::update_shift(std::size_t size)
{
    // New (i.e. shifted) orders. We will have at most 3*N_states new states
    this->orders.resize(3*size);
    this->locations.clear();
    this->locations.reserve(3*size);
    
    // Make sure k=0 is in the first position.
    this->orders[0] = 0;
    this->locations[0] = 0;
    
    // Same for F states.
    this->F.resize(3*size);
    // NOTE memset is slightly faster than std::fill, but requires that Complex
    // has a trivial representation.
    static_assert(
        std::is_trivially_copyable<Complex>::value, 
        "Complex cannot be used with memset");
    std::memset(this->F.data(), 0, this->F.size()*sizeof(decltype(this->F[0])));
    
    // Same for F* states.
    this->F_star.resize(3*size);
    std::memset(
        this->F_star.data(), 0, 
        this->F_star.size()*sizeof(decltype(this->F_star[0])));
    
    // Same for Z states.
    this->Z.resize(3*size);
    std::memset(this->Z.data(), 0, this->Z.size()*sizeof(decltype(this->Z[0])));
}

void
Discrete::Cache
::update_diffusion(
    std::size_t size, decltype(Discrete::_orders) const & orders, 
    Real bin_width)
{
    this->k.resize(size);
    for(std::size_t order=0; order != size; ++order)
    {
        this->k[order] = orders[order]*bin_width;
    }
}

std::size_t
Discrete::Cache
::get_location(long long order) 
{
    auto const location = this->locations.size();
    auto const insert_result = this->locations.try_emplace(order, location);
    if(insert_result.second)
    {
        this->orders[location] = order;
    }
    return insert_result.first->second;
}

}

}
