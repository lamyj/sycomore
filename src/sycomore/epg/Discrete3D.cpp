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
    for(; it!=this->_orders.end(); it+=3)
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
    for(unsigned int order=0; order<this->size(); ++order)
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
        for(std::size_t source=1; source<this->size(); ++source)
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

        this->_orders.resize(3*(destination+1));
        this->_F.resize(destination+1);
        this->_F_star.resize(destination+1);
        this->_Z.resize(destination+1);
        
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

/**
 * @brief Return the location of given order in the states vectors, create it 
 * if missing.
 */
std::size_t get_location(
    std::map<std::array<int64_t, 3>, std::size_t> & locations,
    std::vector<int64_t, xsimd::aligned_allocator<int64_t, 64>> & orders,
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> & F,
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> & F_star,
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> & Z,
    std::array<int64_t, 3> const & order) 
{
    auto iterator = locations.find(order);
    if(iterator == locations.end())
    {
        F.push_back(0.);
        F_star.push_back(0.);
        Z.push_back(0.);
        std::copy(order.begin(), order.end(), std::back_inserter(orders));
        iterator = locations.emplace(order, F.size()-1).first;
    }
    return iterator->second;
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
    
    // New (i.e. shifted) orders. Make sure k=0 is in the first position.
    decltype(this->_orders) orders; 
    orders.reserve(this->_orders.size());
    orders.push_back(0); orders.push_back(0); orders.push_back(0);
    // Same for F states.
    decltype(this->_F) F;
    F.reserve(this->_F.size());
    F.push_back(0.);
    // Same for F* states.
    decltype(this->_F_star) F_star;
    F_star.reserve(this->_F_star.size());
    F_star.push_back(0.);
    // Same for Z states.
    decltype(this->_Z) Z; 
    Z.reserve(this->_Z.size());
    Z.push_back(0.);
    
    // Mapping between a normalized (i.e. folded) order and its location in the
    // states vectors.
    // NOTE: unordered_map is slower than map here, even with a simple hash.
    std::map<Bin, std::size_t> locations;
    locations[{0,0,0}] = 0;
    
    for(std::size_t i=0; i<this->size(); ++i)
    {
        Bin const k{
            this->_orders[3*i+0], this->_orders[3*i+1], this->_orders[3*i+2]};
        if(this->_Z[i] != 0.)
        {
            auto const location = get_location(
                locations, orders, F, F_star, Z, k);
            Z[location] = this->_Z[i];
        }
        
        if(this->_F[i] != 0.)
        {
            Bin k_F{k[0]+delta_k[0], k[1]+delta_k[1], k[2]+delta_k[2]};
            
            // Depending on whether the new order changed half space, conjugate
            // the state and store it in F* instead of F.
            auto destination = &F;
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
                destination = &F_star;
                value = std::conj(value);
            }
            auto const location = get_location(
                locations, orders, F, F_star, Z, k_F);
            (*destination)[location] = value;
        }
        
        // WARNING: F* state at echo is a duplicate of F state.
        if(i != 0 && this->_F_star[i] != 0.)
        {
            // The F* order corresponding to F order k+Δk is -(-k+Δk), i.e. k-Δk
            Bin k_F_star{k[0]-delta_k[0], k[1]-delta_k[1], k[2]-delta_k[2]};
            
            // Same as above.
            auto destination = &F_star;
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
                destination = &F;
                value = std::conj(value);
            }
            auto const location = get_location(
                locations, orders, F, F_star, Z, k_F_star);
            (*destination)[location] = value;
        }
    }
    
    // Update the current orders and states with the new ones.
    this->_orders = std::move(orders);
    this->_F = std::move(F);
    this->_F_star = std::move(F_star);
    this->_Z = std::move(Z);
    
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
    bool all_zero=true;
    for(std::size_t i=0; i<9; ++i)
    {
        if(this->species.get_D()[i].magnitude != 0)
        {
            all_zero = false;
            break;
        }
    }
    if(all_zero)
    {
        return;
    }
    
    auto const tau = duration.magnitude;
    auto const bin_width = this->_bin_width.magnitude;
    
    using AlignedVector = std::vector<Real, xsimd::aligned_allocator<Real, 64>>;
    
    std::vector<AlignedVector> k(3);
    k[0].resize(this->size());
    k[1].resize(this->size());
    k[2].resize(this->size());
    for(std::size_t order=0; order<this->size(); ++order)
    {
        k[0][order] = this->_orders[3*order+0]*bin_width;
        k[1][order] = this->_orders[3*order+1]*bin_width;
        k[2][order] = this->_orders[3*order+2]*bin_width;
    }
    
    std::vector<Real, xsimd::aligned_allocator<Real, 64>> const delta_k{
        sycomore::gamma.magnitude*gradient[0].magnitude*tau,
        sycomore::gamma.magnitude*gradient[1].magnitude*tau,
        sycomore::gamma.magnitude*gradient[2].magnitude*tau
    };
    
    AlignedVector b_L_D(this->size(), 0.);
    AlignedVector b_T_plus_D(this->size(), 0.);
    AlignedVector b_T_minus_D(this->size(), 0.);
    for(std::size_t m=0; m<3; ++m)
    {
        for(std::size_t n=0; n<3; ++n)
        {
            auto const delta_k_product_term = 
                1./3. * tau * delta_k[m] * delta_k[n];
            
            simd_api::diffusion_3d_b(
                k[m].data(), k[n].data(), delta_k[m], delta_k[n], 
                delta_k_product_term, tau, 
                this->species.get_D()[3*m+n].magnitude, 
                b_L_D.data(), b_T_plus_D.data(), b_T_minus_D.data(),
                this->_F.size());
        }
    }
    
    simd_api::diffusion_3d(
        b_L_D.data(), b_T_plus_D.data(), b_T_minus_D.data(),
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

}

}
