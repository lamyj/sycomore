#include "Discrete3D.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <sstream>
#include <utility>
#include <vector>

#include <xtensor/xio.hpp>

#include "sycomore/Array.h"
#include "sycomore/epg/Base.h"
#include "sycomore/epg/robin_hood.h"
#include "sycomore/epg/operators.h"
#include "sycomore/epg/simd_api.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/sycomore.h"

namespace sycomore
{

namespace epg
{

Discrete3D
::Discrete3D(
    Species const & species, Vector3R const & initial_magnetization,
    Quantity bin_width)
: Base(species, initial_magnetization, 1),
    _orders{{0,0,0}}, _bin_width(bin_width), _cache(this->_model.pools)
{
    // Nothing else.
}

Discrete3D
::Discrete3D(
    Species const & species_a, Species const & species_b,
    Vector3R const & M0_a, Vector3R const & M0_b,
    Quantity const & k_a, Quantity const & delta_b, Quantity bin_width)
: Base(species_a, species_b, M0_a, M0_b, k_a, delta_b, 1),
    _orders{{0,0,0}}, _bin_width(bin_width), _cache(this->_model.pools)
{
    // Nothing else
}

Discrete3D
::Discrete3D(
    Species const & species_a, Quantity const & R1_b_or_T1_b,
    Vector3R const & M0_a, Vector3R const & M0_b,
    Quantity const & k_a, Quantity bin_width)
: Base(species_a, R1_b_or_T1_b, M0_a, M0_b, k_a, 1),
    _orders{{0,0,0}}, _bin_width(bin_width), _cache(this->_model.pools)
{
    // Nothing else
}

std::size_t
Discrete3D
::size() const
{
    return this->_orders.size()/3;
}

TensorQ<2>
Discrete3D
::orders() const
{
    TensorQ<2> orders(TensorQ<2>::shape_type{this->size(), 3});
    std::transform(
        this->_orders.begin(), this->_orders.end(), orders.begin(),
        [&](Orders::value_type const & k){ return k*this->_bin_width; });
    return orders;
}

ArrayC
Discrete3D
::state(Order const & order) const
{
    if(order.size() != 3)
    {
        std::ostringstream message;
        message << "Order must have 3 elements, not " << order.size();
        throw std::runtime_error(message.str());
    }
    
    Bin bin{
        static_cast<int64_t>(std::round(order[0]/this->_bin_width)),
        static_cast<int64_t>(std::round(order[1]/this->_bin_width)),
        static_cast<int64_t>(std::round(order[2]/this->_bin_width)) };
    auto it = this->_orders.begin();
    for(auto end = this->_orders.end(); it != end; it+=3)
    {
        if(bin == Bin{it[0], it[1], it[2]})
        {
            break;
        }
    }
    if(it == this->_orders.end())
    {
        std::ostringstream message;
        message << "No such order: " << order;
        throw std::runtime_error(message.str());
    }

    std::size_t const position = (it-this->_orders.begin())/3;
    return this->state(position);
}

void
Discrete3D
::apply_time_interval(
    Quantity const & duration, Vector3Q const & gradient)
{
    if(duration.magnitude == 0)
    {
        return;
    }

    this->relaxation(duration);
    this->diffusion(duration, gradient);
    this->shift(duration, gradient);
    this->off_resonance(duration);
    
    this->_elapsed += duration.magnitude;
    
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
Discrete3D
::apply_time_interval(TimeInterval const & interval)
{
    this->apply_time_interval(
        interval.get_duration(),
        {
            interval.get_gradient_amplitude()[0],
            interval.get_gradient_amplitude()[1],
            interval.get_gradient_amplitude()[2]
        });
}

void
Discrete3D
::shift(Quantity const & duration, Vector3Q const & gradient)
{
    // Compute dephasing and return early if it is null.
    Bin const delta_k {
        static_cast<int64_t>(std::round(
            sycomore::gamma.magnitude*gradient[0].magnitude*duration.magnitude
            / this->_bin_width.magnitude)),
        static_cast<int64_t>(std::round(
            sycomore::gamma.magnitude*gradient[1].magnitude*duration.magnitude
            / this->_bin_width.magnitude)),
        static_cast<int64_t>(std::round(
            sycomore::gamma.magnitude*gradient[2].magnitude*duration.magnitude
            / this->_bin_width.magnitude))};
    if(delta_k[0] == 0 && delta_k[1] == 0 && delta_k[2] == 0)
    {
        return;
    }
    
    this->_cache.update_shift(this->size());
    
    for(std::size_t i=0, end=this->size(); i != end; ++i)
    {
        Bin const k{
            this->_orders[3*i+0], this->_orders[3*i+1], this->_orders[3*i+2]};
        
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
                Bin k_F{k[0]+delta_k[0], k[1]+delta_k[1], k[2]+delta_k[2]};
                
                // Depending on whether the new order changed half space, conjugate
                // the state and store it in F* instead of F.
                auto destination = &this->_cache.F[pool];
                auto & value = this->_model.F[pool][i];
                // WARNING: the half-space (k_F[0] >= 0) contains conjugate
                // states, e.g. ([0, y, 0], [0, -y, 0]) or more generally all
                // pairs of the form ([0, y, z], [0, -y, -z]). Solve this by
                // including only part of the y and z axes.
                if(
                    k_F[0] < 0 
                    || (k_F[0] == 0 && k_F[1] < 0)
                    || (k_F[0] == 0 && k_F[1] == 0 && k_F[2] <0))
                {
                    k_F[0] *= -1;
                    k_F[1] *= -1;
                    k_F[2] *= -1;
                    destination = &this->_cache.F_star[pool];
                    value = std::conj(value);
                }
                (*destination)[this->_cache.get_location(k_F)] = value;
            }
            
            // WARNING: F* state at echo is a duplicate of F state.
            state = this->_model.F_star[pool][i];
            if(i != 0 && state != 0.)
            {
                // The F* order corresponding to F order k+Δk is -(-k+Δk),
                // i.e. k-Δk
                Bin k_F_star{k[0]-delta_k[0], k[1]-delta_k[1], k[2]-delta_k[2]};
                
                // Same as above.
                auto destination = &this->_cache.F_star[pool];
                auto & value = this->_model.F_star[pool][i];
                // cf. WARNING about conjugation in F case.
                if(
                    k_F_star[0] < 0 
                    || (k_F_star[0] == 0 && k_F_star[1] < 0)
                    || (k_F_star[0] == 0 && k_F_star[1] == 0 && k_F_star[2] <0))
                {
                    k_F_star[0] *= -1;
                    k_F_star[1] *= -1;
                    k_F_star[2] *= -1;
                    destination = &this->_cache.F[pool];
                    value = std::conj(value);
                }
                (*destination)[this->_cache.get_location(k_F_star)] = value;
            }
        }
    }
    
    // Update the current orders and states with the new ones.
    this->_cache.orders.resize(this->_cache.locations.size()*3);
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
Discrete3D
::diffusion(Quantity const & duration, Vector3Q const & gradient)
{
    if(
        std::all_of(
            this->_model.species.begin(), this->_model.species.end(),
            [](Species const & s) {
                return std::all_of(
                    s.get_D().begin(), s.get_D().end(),
                    [](Quantity const & x) { return x.magnitude == 0;});
            }
        ))
    {
        return;
    }
    
    auto const tau = duration.magnitude;
    
    Cache::RealVector const delta_k{
        sycomore::gamma.magnitude*gradient[0].magnitude*tau,
        sycomore::gamma.magnitude*gradient[1].magnitude*tau,
        sycomore::gamma.magnitude*gradient[2].magnitude*tau
    };
    if(std::all_of(delta_k.begin(), delta_k.end(), [](Real x) { return x == 0;}))
    {
        return;
    }
    
    this->_cache.update_diffusion(
        this->size(), this->_orders, this->_bin_width.magnitude);
    
    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
    {
        auto & species = this->_model.species[pool];
        auto & F = this->_model.F[pool];
        auto & F_star = this->_model.F_star[pool];
        auto & Z = this->_model.Z[pool];
        
        auto & cache = this->_cache;
        
        for(std::size_t m=0; m<3; ++m)
        {
            for(std::size_t n=0; n<3; ++n)
            {
                auto const delta_k_product_term = 
                    1./3. * tau * delta_k[m] * delta_k[n];
                
                simd_api::diffusion_3d_b(
                    cache.k[m].data(), cache.k[n].data(), 
                    delta_k[m], delta_k[n], delta_k_product_term,
                    tau, species.get_D().unchecked(m, n).magnitude,
                    cache.b_L_D.data(), 
                    cache.b_T_plus_D.data(), cache.b_T_minus_D.data(),
                    F.size());
            }
        }
        
        simd_api::diffusion_3d(
            cache.b_L_D.data(), 
            cache.b_T_plus_D.data(), cache.b_T_minus_D.data(),
            F.data(), F_star.data(), Z.data(), this->size());
    }
    
}

Quantity const & 
Discrete3D
::bin_width() const
{
    return this->_bin_width;
}

Discrete3D::Cache
::Cache(std::size_t pools)
: orders(0), F(pools), F_star(pools), Z(pools)
{
    // Nothing else.
}

void
Discrete3D::Cache
::update_shift(std::size_t size)
{
    // New (i.e. shifted) orders. We will have at most 3*N_states new states
    this->orders.resize(3*3*size);
    this->locations.clear();
    this->locations.reserve(3*size);
    
    // Make sure k=0 is in the first position.
    this->orders[0] = this->orders[1] = this->orders[2] = 0;
    this->locations[{0,0,0}] = 0;
    
    // Same for F states.
    for(auto & F: this->F)
    {
        F.resize(3*size);
        // NOTE memset is slightly faster than std::fill, but requires that Complex
        // has a trivial representation.
        static_assert(
            std::is_trivially_copyable<Complex>::value, 
            "Complex cannot be used with memset");
        std::memset(
            reinterpret_cast<void*>(F.data()), 0, F.size()*sizeof(Complex));
    }
    
    // Same for F* states.
    for(auto & F_star: this->F_star)
    {
        F_star.resize(3*size);
        std::memset(
            reinterpret_cast<void*>(F_star.data()), 0,
            F_star.size()*sizeof(Complex));
    }
    
    // Same for Z states.
    for(auto & Z: this->Z)
    {
        Z.resize(3*size);
        std::memset(
            reinterpret_cast<void*>(Z.data()), 0, Z.size()*sizeof(Complex));
    }
}

void
Discrete3D::Cache
::update_diffusion(std::size_t size, Orders const & orders, Real bin_width)
{
    this->k[0].resize(size);
    this->k[1].resize(size);
    this->k[2].resize(size);
    for(std::size_t order=0; order != size; ++order)
    {
        this->k[0][order] = orders[3*order+0]*bin_width;
        this->k[1][order] = orders[3*order+1]*bin_width;
        this->k[2][order] = orders[3*order+2]*bin_width;
    }
    
    this->b_L_D.resize(size);
    std::fill(this->b_L_D.begin(), this->b_L_D.end(), 0);
    
    this->b_T_plus_D.resize(size);
    std::fill(this->b_T_plus_D.begin(), this->b_T_plus_D.end(), 0);
    
    this->b_T_minus_D.resize(size);
    std::fill(this->b_T_minus_D.begin(), this->b_T_minus_D.end(), 0);
}

std::size_t
Discrete3D::Cache
::get_location(Bin const & order) 
{
    auto const location = this->locations.size();
    auto const insert_result = this->locations.try_emplace(order, location);
    if(insert_result.second)
    {
        std::copy(order.begin(), order.end(), this->orders.data()+3*location);
    }
    return insert_result.first->second;
}

}

}
