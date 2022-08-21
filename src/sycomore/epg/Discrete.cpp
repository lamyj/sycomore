#include "Discrete.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/epg/Base.h"
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
    Quantity bin_width)
: Base(species, initial_magnetization, 1),
    _bin_width(bin_width), _orders{0}, _cache(this->_model.pools)
{
    // Nothing else
}

Discrete
::Discrete(
    Species const & species_a, Species const & species_b,
    Magnetization const & M0_a, Magnetization const & M0_b,
    Quantity const & k_a, Quantity const & delta_b, Quantity bin_width)
: Base(species_a, species_b, M0_a, M0_b, k_a, delta_b, 1),
    _bin_width(bin_width), _orders{0}, _cache(this->_model.pools)
{
    // Nothing else
}

Discrete
::Discrete(
    Species const & species_a, Quantity const & R1_b_or_T1_b,
    Magnetization const & M0_a, Magnetization const & M0_b,
    Quantity const & k_a, Quantity bin_width)
: Base(species_a, R1_b_or_T1_b, M0_a, M0_b, k_a, 1),
    _bin_width(bin_width), _orders{0}, _cache(this->_model.pools)
{
    // Nothing else
}

std::size_t
Discrete
::size() const
{
    return this->_orders.size();
}

std::vector<Discrete::Order>
Discrete
::orders() const
{
    std::vector<Order> orders(this->size());
    std::transform(
        this->_orders.begin(), this->_orders.end(), orders.begin(),
        [&](Orders::value_type const & k){ return k*this->_bin_width; });
    return orders;
}

std::vector<Complex>
Discrete
::state(Order const & order) const
{
    std::size_t const k = std::round(double(order/this->_bin_width));

    auto const it = std::find(this->_orders.begin(), this->_orders.end(), k);
    if(it == this->_orders.end())
    {
        std::ostringstream message;
        message << "No such order: " << order;
        throw std::runtime_error(message.str());
    }

    std::size_t const position = it-this->_orders.begin();
    return this->state(position);
}

void
Discrete
::apply_time_interval(
    Quantity const & duration, Quantity const & gradient)
{
    if(duration.magnitude == 0)
    {
        return;
    }
    
    this->relaxation(duration);
    this->diffusion(duration, gradient);
    this->bulk_motion(duration, gradient);
    this->shift(duration, gradient);
    this->off_resonance(duration);
    
    if(this->threshold > 0)
    {
        auto const threshold_squared = std::pow(this->threshold, 2);
        
        // Always include the zero order (implicit since we start at 1),
        // include other order if population is above threshold.
        std::size_t destination=1;
        for(std::size_t source=1, end=this->size(); source != end; ++source)
        {
            Real max_magnitude_squared = 0.;
            for(std::size_t pool=0; pool<this->_model.pools; ++pool)
            {
                using std::pow; using std::abs;
                auto const magnitude_squared = 
                    pow(abs(this->_model.F[pool][source]), 2)
                    +pow(abs(this->_model.F_star[pool][source]), 2)
                    +pow(abs(this->_model.Z[pool][source]), 2);
                max_magnitude_squared = std::max(
                    max_magnitude_squared, magnitude_squared);
            }
            
            if(max_magnitude_squared >= threshold_squared)
            {
                if(source != destination)
                {
                    this->_orders[destination] = this->_orders[source];
                    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
                    {
                        this->_model.F[pool][destination] =
                            this->_model.F[pool][source];
                        this->_model.F_star[pool][destination] =
                            this->_model.F_star[pool][source];
                        this->_model.Z[pool][destination] =
                            this->_model.Z[pool][source];
                    }
                }
                ++destination;
            }
        }
        
        this->_orders.resize(destination);
        for(std::size_t pool=0; pool<this->_model.pools; ++pool)
        {
            this->_model.F[pool].resize(destination);
            this->_model.F_star[pool].resize(destination);
            this->_model.Z[pool].resize(destination);
                
            // No need to update the iterator pointing to the echo magnetization.
        }
    }
}

void
Discrete
::apply_time_interval(TimeInterval const & interval)
{
    this->apply_time_interval(
        interval.get_duration(), interval.get_gradient_amplitude()[0]);
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
        
        for(std::size_t pool=0; pool<this->_model.pools; ++pool)
        {
            auto state = this->_model.Z[pool][i];
            if(state != 0.)
            {
                this->_cache.Z[pool][this->_cache.get_location(k)] = state;
            }
            
            state = this->_model.F[pool][i];
            if(state != 0.)
            {
                auto k_F = k+delta_k;
                
                // Depending on whether the new order changed half space,
                // conjugate the state and store it in F* instead of F.
                auto destination = &this->_cache.F[pool];
                if(k_F < 0)
                {
                    k_F *= -1;
                    destination = &this->_cache.F_star[pool];
                    state = std::conj(state);
                }
                (*destination)[this->_cache.get_location(k_F)] = state;
            }
            
            // WARNING: F* state at echo is a duplicate of F state.
            state = this->_model.F_star[pool][i];
            if(i != 0 && state != 0.)
            {
                // The F* order corresponding to F order k+Δk is -(-k+Δk),
                // i.e. k-Δk
                auto k_F_star = k-delta_k;
                
                // Same as above.
                auto destination = &this->_cache.F_star[pool];
                if(k_F_star < 0)
                {
                    k_F_star *= -1;
                    destination = &this->_cache.F[pool];
                    state = std::conj(state);
                }
                (*destination)[this->_cache.get_location(k_F_star)] = state;
            }
        }
    }
    
    // Update the current orders and states with the new ones.
    this->_cache.orders.resize(this->_cache.locations.size());
    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
    {
        this->_cache.F[pool].resize(this->_cache.locations.size());
        this->_cache.F_star[pool].resize(this->_cache.locations.size());
        this->_cache.Z[pool].resize(this->_cache.locations.size());
    }
        
    // Use swap and not move since we keep the temporary variables between
    // runs.
    std::swap(this->_cache.orders, this->_orders);
    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
    {
        std::swap(this->_cache.F[pool], this->_model.F[pool]);
        std::swap(this->_cache.F_star[pool], this->_model.F_star[pool]);
        std::swap(this->_cache.Z[pool], this->_model.Z[pool]);
    }
    
    // Update the conjugate states of the echo magnetization.
    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
    {
        if(this->_model.F[pool][0] != 0.)
        {
            this->_model.F_star[pool][0] = std::conj(this->_model.F[pool][0]);
        }
        else if(this->_model.F_star[pool][0] != 0.)
        {
            this->_model.F[pool][0] = std::conj(this->_model.F_star[pool][0]);
        }
    }
}

void
Discrete
::diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(
        std::all_of(
            this->_model.species.begin(), this->_model.species.end(),
            [](Species const & s) { return s.get_D()[0].magnitude == 0; }))
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
    
    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
    {
        auto const & D = this->_model.species[pool].get_D()[0].magnitude;
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
        this->_model, this->size());
}

Quantity const & 
Discrete
::bin_width() const
{
    return this->_bin_width;
}

Discrete::Cache
::Cache(std::size_t pools)
: orders(0), F(pools), F_star(pools), Z(pools)
{
    // Nothing else.
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
    for(auto & F: this->F)
    {
        F.resize(3*size);
        // NOTE memset is slightly faster than std::fill, but requires that
        // Complex has a trivial representation.
        static_assert(
            std::is_trivially_copyable<Complex>::value, 
            "Complex cannot be used with memset");
        std::memset(F.data(), 0, F.size()*sizeof(Complex));
    }
    
    // Same for F* states.
    for(auto & F_star: this->F_star)
    {
        F_star.resize(3*size);
        std::memset(F_star.data(), 0, F_star.size()*sizeof(Complex));
    }
    
    // Same for Z states.
    for(auto & Z: this->Z)
    {
        Z.resize(3*size);
        std::memset(Z.data(), 0, Z.size()*sizeof(Complex));
    }
}

void
Discrete::Cache
::update_diffusion(std::size_t size, Orders const & orders, Real bin_width)
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
