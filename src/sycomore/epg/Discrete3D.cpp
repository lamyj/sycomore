#include "Discrete3D.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include "sycomore/Array.h"
#include "sycomore/epg/Base.h"
#include "sycomore/epg/pool_model.h"
#include "sycomore/epg/pool_storage.h"
#include "sycomore/epg/robin_hood.h"
#include "sycomore/epg/operators.h"
#include "sycomore/epg/simd_api.h"
#include "sycomore/magnetization.h"
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
    Species const & species, Magnetization const & initial_magnetization,
    Quantity bin_width, Real threshold)
: Base(species, initial_magnetization, 1), _bin_width(bin_width)
{
    this->threshold = threshold;
    this->_orders = {0,0,0};
}

std::size_t
Discrete3D
::size() const
{
    return this->_storage->F.size();
}

std::vector<Quantity>
Discrete3D
::orders() const
{
    std::vector<Quantity> orders(this->_orders.size());
    std::transform(
        this->_orders.begin(), this->_orders.end(), orders.begin(),
        [&](decltype(this->_orders[0]) k){ return k*this->_bin_width; });
    return orders;
}

std::vector<Complex>
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
        static_cast<int64_t>(order[0]/this->_bin_width),
        static_cast<int64_t>(order[1]/this->_bin_width),
        static_cast<int64_t>(order[2]/this->_bin_width) };
    auto it=this->_orders.begin();
    for(auto end=this->_orders.end(); it != end; it+=3)
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

    auto const index = (it-this->_orders.begin())/3;
    return this->Base::state(index);
}

void
Discrete3D
::apply_time_interval(
    Quantity const & duration, Array<Quantity> const & gradient)
{
    if(duration.magnitude == 0)
    {
        return;
    }

    this->relaxation(duration);
    this->diffusion(duration, gradient);
    this->shift(duration, gradient);
    this->off_resonance(duration);
    
    if(this->threshold > 0)
    {
        auto const threshold_squared = std::pow(this->threshold, 2);
        
        std::size_t destination=1;
        for(std::size_t source=1, end=this->size(); source != end; ++source)
        {
            auto const magnitude_squared =
                std::pow(std::abs(this->_storage->F[source]), 2)
                +std::pow(std::abs(this->_storage->F_star[source]), 2)
                +std::pow(std::abs(this->_storage->Z[source]), 2);
            // Always include the zero order (implicit since we start at 1),
            // include other order if population is above threshold.
            if(magnitude_squared >= threshold_squared)
            {
                if(source != destination)
                {
                    std::copy(
                        this->_orders.data()+3*source, 
                        this->_orders.data()+3*(1+source), 
                        this->_orders.data()+3*destination);
                    this->_storage->F[destination] = this->_storage->F[source];
                    this->_storage->F_star[destination] = this->_storage->F_star[source];
                    this->_storage->Z[destination] = this->_storage->Z[source];
                }
                ++destination;
            }
        }

        this->_orders.resize(3*(destination));
        this->_storage->F.resize(destination);
        this->_storage->F_star.resize(destination);
        this->_storage->Z.resize(destination);
        
        // No need to update the iterator pointing to the echo magnetization.
    }
}

void
Discrete3D
::apply_time_interval(TimeInterval const & interval)
{
    this->apply_time_interval(
        interval.get_duration(), interval.get_gradient_amplitude());
}

void
Discrete3D
::shift(Quantity const & duration, Array<Quantity> const & gradient)
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
        if(this->_storage->Z[i] != 0.)
        {
            this->_cache.Z[this->_cache.get_location(k)] = this->_storage->Z[i];
        }
        
        if(this->_storage->F[i] != 0.)
        {
            Bin k_F{k[0]+delta_k[0], k[1]+delta_k[1], k[2]+delta_k[2]};
            
            // Depending on whether the new order changed half space, conjugate
            // the state and store it in F* instead of F.
            auto destination = &this->_cache.F;
            auto & value = this->_storage->F[i];
            // WARNING: the half-space (k_F[0] >= 0) contains conjugate states,
            // e.g. ([0, y, 0], [0, -y, 0]) or more generally all pairs of the
            // form ([0, y, z], [0, -y, -z]). Solve this by including only part
            // of the y and z axes.
            if(
                k_F[0] < 0 
                || (k_F[0] == 0 && k_F[1] < 0)
                || (k_F[0] == 0 && k_F[1] == 0 && k_F[2] <0))
            {
                k_F[0] *= -1;
                k_F[1] *= -1;
                k_F[2] *= -1;
                destination = &this->_cache.F_star;
                value = std::conj(value);
            }
            (*destination)[this->_cache.get_location(k_F)] = value;
        }
        
        // WARNING: F* state at echo is a duplicate of F state.
        if(i != 0 && this->_storage->F_star[i] != 0.)
        {
            // The F* order corresponding to F order k+Δk is -(-k+Δk), i.e. k-Δk
            Bin k_F_star{k[0]-delta_k[0], k[1]-delta_k[1], k[2]-delta_k[2]};
            
            // Same as above.
            auto destination = &this->_cache.F_star;
            auto & value = this->_storage->F_star[i];
            // cf. WARNING about conjugation in F case.
            if(
                k_F_star[0] < 0 
                || (k_F_star[0] == 0 && k_F_star[1] < 0)
                || (k_F_star[0] == 0 && k_F_star[1] == 0 && k_F_star[2] <0))
            {
                k_F_star[0] *= -1;
                k_F_star[1] *= -1;
                k_F_star[2] *= -1;
                destination = &this->_cache.F;
                value = std::conj(value);
            }
            (*destination)[this->_cache.get_location(k_F_star)] = value;
        }
    }
    
    // Update the current orders and states with the new ones.
    this->_cache.orders.resize(this->_cache.locations.size()*3);
    this->_cache.F.resize(this->_cache.locations.size());
    this->_cache.F_star.resize(this->_cache.locations.size());
    this->_cache.Z.resize(this->_cache.locations.size());
    
    // Use swap and not move since we keep the temporary variables between runs.
    std::swap(this->_cache.orders, this->_orders);
    std::swap(this->_cache.F, this->_storage->F);
    std::swap(this->_cache.F_star, this->_storage->F_star);
    std::swap(this->_cache.Z, this->_storage->Z);
    
    // Update the conjugate states of the echo magnetization.
    if(this->_storage->F[0] != 0.)
    {
        this->_storage->F_star[0] = std::conj(this->_storage->F[0]);
    }
    else if(this->_storage->F_star[0] != 0.)
    {
        this->_storage->F[0] = std::conj(this->_storage->F_star[0]);
    }
}

void
Discrete3D
::diffusion(Quantity const & duration, Array<Quantity> const & gradient)
{
    if(std::all_of(
        this->_model->species.get_D().begin(),
        this->_model->species.get_D().end(), 
        [](Quantity const & x) { return x.magnitude == 0;}))
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
    
    for(std::size_t m=0; m<3; ++m)
    {
        for(std::size_t n=0; n<3; ++n)
        {
            auto const delta_k_product_term = 
                1./3. * tau * delta_k[m] * delta_k[n];
            
            simd_api::diffusion_3d_b(
                this->_cache.k[m].data(), this->_cache.k[n].data(), 
                delta_k[m], delta_k[n], delta_k_product_term,
                tau, this->_model->species.get_D()[3*m+n].magnitude,
                this->_cache.b_L_D.data(), 
                this->_cache.b_T_plus_D.data(), this->_cache.b_T_minus_D.data(),
                this->_storage->F.size());
        }
    }
    
    simd_api::diffusion_3d(
        this->_cache.b_L_D.data(), 
        this->_cache.b_T_plus_D.data(), this->_cache.b_T_minus_D.data(),
        this->_storage->F.data(),
        this->_storage->F_star.data(),
        this->_storage->Z.data(), 
        this->_storage->F.size());
}

Quantity const & 
Discrete3D
::bin_width() const
{
    return this->_bin_width;
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
Discrete3D::Cache
::update_diffusion(
    std::size_t size, decltype(Discrete3D::_orders) const & orders, 
    Real bin_width)
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
