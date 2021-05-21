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
: species(species), _bin_width(bin_width), _F(1, 0), _F_star(1, 0), _Z(1, 0), 
    threshold(threshold)
{
    auto const M = as_complex_magnetization(initial_magnetization);
    this->_orders = {0,0,0};
    this->_F[0] = std::sqrt(2)*M.p;
    this->_F_star[0] = std::sqrt(2)*M.m;
    this->_Z[0] = M.z;
    this->_M_z_eq = M.z;
}

std::size_t
Discrete3D
::size() const
{
    return this->_F.size();
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

Discrete3D::State
Discrete3D
::state(std::size_t order) const
{
    return {this->_F[order], this->_F_star[order], this->_Z[order]};
}

Discrete3D::State
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
    return this->state(index);
}

std::vector<Complex>
Discrete3D
::states() const
{
    std::vector<Complex> result(3*this->size());
    for(std::size_t order=0, end=this->size(); order != end; ++order)
    {
        result[3*order+0] = this->_F[order];
        result[3*order+1] = this->_F_star[order];
        result[3*order+2] = this->_Z[order];
    }
    return result;
}

Complex const &
Discrete3D
::echo() const
{
    return this->_F[0];
}

void
Discrete3D
::apply_pulse(Quantity angle, Quantity phase)
{
    simd_api::apply_pulse(
        operators::pulse(angle.magnitude, phase.magnitude), 
        this->_F.data(), this->_F_star.data(), this->_Z.data(), 
        this->_F.size());
}

void
Discrete3D
::apply_time_interval(
    Quantity const & duration, Array<Quantity> const & gradient, Real threshold)
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
                    std::copy(
                        this->_orders.data()+3*source, 
                        this->_orders.data()+3*(1+source), 
                        this->_orders.data()+3*destination);
                    this->_F[destination] = this->_F[source];
                    this->_F_star[destination] = this->_F_star[source];
                    this->_Z[destination] = this->_Z[source];
                }
                ++destination;
            }
        }

        this->_orders.resize(3*(destination));
        this->_F.resize(destination);
        this->_F_star.resize(destination);
        this->_Z.resize(destination);
        
        // No need to update the iterator pointing to the echo magnetization.
    }
}

void
Discrete3D
::apply_time_interval(TimeInterval const & interval)
{
    this->apply_time_interval(
        interval.get_duration(), interval.get_gradient_amplitude(), 0);
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
        if(this->_Z[i] != 0.)
        {
            this->_cache.Z[this->_cache.get_location(k)] = this->_Z[i];
        }
        
        if(this->_F[i] != 0.)
        {
            Bin k_F{k[0]+delta_k[0], k[1]+delta_k[1], k[2]+delta_k[2]};
            
            // Depending on whether the new order changed half space, conjugate
            // the state and store it in F* instead of F.
            auto destination = &this->_cache.F;
            auto & value = this->_F[i];
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
        if(i != 0 && this->_F_star[i] != 0.)
        {
            // The F* order corresponding to F order k+Δk is -(-k+Δk), i.e. k-Δk
            Bin k_F_star{k[0]-delta_k[0], k[1]-delta_k[1], k[2]-delta_k[2]};
            
            // Same as above.
            auto destination = &this->_cache.F_star;
            auto & value = this->_F_star[i];
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
Discrete3D
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
        this->_F.size());
    
    this->_Z[0] += this->_M_z_eq*(1.-E.first);
}

void
Discrete3D
::diffusion(Quantity const & duration, Array<Quantity> const & gradient)
{
    if(std::all_of(
        this->species.get_D().begin(), this->species.get_D().end(), 
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
                tau, this->species.get_D()[3*m+n].magnitude,
                this->_cache.b_L_D.data(), 
                this->_cache.b_T_plus_D.data(), this->_cache.b_T_minus_D.data(),
                this->_F.size());
        }
    }
    
    simd_api::diffusion_3d(
        this->_cache.b_L_D.data(), 
        this->_cache.b_T_plus_D.data(), this->_cache.b_T_minus_D.data(),
        this->_F.data(), this->_F_star.data(), this->_Z.data(), 
        this->_F.size());
}

void
Discrete3D
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
            this->_F.size());
    }
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
