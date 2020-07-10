#include "Discrete.h"

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "sycomore/epg/operators.h"
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
    // Store magnetization as lines of F̃, F̃^*_, Z̃
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_F[0] = std::sqrt(2)*magnetization.p;
    this->_F_star[0] = std::sqrt(2)*magnetization.m;
    this->_Z[0] = magnetization.z;
    
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
    auto const T = operators::pulse(angle.magnitude, phase.magnitude);
    
    #pragma omp parallel for
    for(int order=0; order<this->_orders.size(); ++order)
    {
        auto const & F = this->_F[order];
        auto const & F_star = this->_F_star[order];
        auto const & Z = this->_Z[order];
    
        auto const F_new =      T[3*0+0] * F + T[3*0+1] * F_star + T[3*0+2] * Z;
        auto const F_star_new = T[3*1+0] * F + T[3*1+1] * F_star + T[3*1+2] * Z;
        auto const Z_new =      T[3*2+0] * F + T[3*2+1] * F_star + T[3*2+2] * Z;
    
        this->_F[order] = F_new;
        this->_F_star[order] = F_star_new;
        this->_Z[order] = Z_new;
    }
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
    this->shift(duration, gradient);
    this->off_resonance(duration);
    
    if(threshold == 0)
    {
        threshold = this->threshold;
    }
    if(threshold > 0)
    {
        auto const threshold_squared = std::pow(threshold, 2);
        decltype(this->_orders) orders; orders.reserve(this->_orders.size());
        decltype(this->_F) F; F.reserve(this->_F.size());
        decltype(this->_F_star) F_star; F_star.reserve(this->_F_star.size());
        decltype(this->_Z) Z; Z.reserve(this->_Z.size());

        for(std::size_t index=0; index<this->_orders.size(); ++index)
        {
            auto const magnitude_squared =
                std::pow(std::abs(this->_F[index]), 2)
                +std::pow(std::abs(this->_F_star[index]), 2)
                +std::pow(std::abs(this->_Z[index]), 2);
            if(magnitude_squared >= threshold_squared)
            {
                orders.push_back(this->_orders[index]);
                F.push_back(this->_F[index]);
                F_star.push_back(this->_F_star[index]);
                Z.push_back(this->_Z[index]);
            }
        }

        this->_orders = std::move(orders);
        this->_F = std::move(F);
        this->_F_star = std::move(F_star);
        this->_Z = std::move(Z);
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
    auto const delta_k = (long long)(std::round(
        double(sycomore::gamma*gradient*duration / this->_bin_width)));
    
    if(delta_k == 0)
    {
        return;
    }
    
    // Unfold the F̃-states
    decltype(this->_orders) F_orders(2*this->_orders.size()-1);
    decltype(this->_F) F(2*this->_orders.size()-1);
    // F̃^+ on the right side
    for(long long i=0; i<this->_orders.size(); ++i)
    {
        F_orders[this->_orders.size()-1 + i] = this->_orders[i];
        F[this->_orders.size()-1 + i] = this->_F[i];
    }
    // F̃^{-*} on the left side, reversed order
    for(long long i=0; i<this->_orders.size(); ++i)
    {
        F_orders[this->_orders.size()-1 - i] = -this->_orders[i];
        F[this->_orders.size()-1 - i] = std::conj(this->_F_star[i]);
    }
    
    // Shift according to Δk
    for(auto && k: F_orders)
    {
        k += delta_k;
    }

    // Z̃-orders do not change
    
    // Fold the F̃-states: build the new orders (using a vector provides quicker
    // iteration later than using a set) …
    decltype(this->_orders) orders(F_orders.size()+this->_orders.size());
    std::transform(
        F_orders.begin(), F_orders.end(), orders.begin(), std::llabs);
    std::copy(
        this->_orders.begin(), this->_orders.end(), 
        orders.begin()+F_orders.size());
    std::sort(orders.begin(), orders.end());
    auto const last = std::unique(orders.begin(), orders.end());
    orders.resize(last-orders.begin());

    // … then build the new states (create a new array for Z̃ states as we are
    // getting their values from this->_Z)…
    this->_F.resize(orders.size()); 
    std::fill(this->_F.begin(), this->_F.end(), 0);
    this->_F_star.resize(orders.size()); 
    std::fill(this->_F_star.begin(), this->_F_star.end(), 0);
    decltype(this->_Z) Z_new(orders.size(), 0);

    // Find the smallest positive F̃-order, move towards larger orders
    auto F_plus_order_it = std::lower_bound(F_orders.cbegin(), F_orders.cend(), 0);
    auto F_plus_state_it = F.cbegin()+(F_plus_order_it-F_orders.cbegin());

    // Find the largest negative F̃-order, move towards smaller orders
    decltype(F_orders)::const_reverse_iterator F_minus_order_it(F_plus_order_it);
    auto F_minus_state_it = F.crbegin()+(F_minus_order_it-F_orders.crbegin());

    auto Z_order_it = this->_orders.cbegin();
    auto Z_state_it = this->_Z.cbegin();

    for(std::size_t index=0; index < orders.size(); ++index)
    {
        // If one of the iterators matches the order, fill the corresponding
        // state and advance

        if(F_plus_order_it != F_orders.cend() && *F_plus_order_it == orders[index])
        {
            this->_F[index] = *F_plus_state_it;
            ++F_plus_order_it;
            ++F_plus_state_it;
        }
        if(F_minus_order_it != F_orders.crend() && *F_minus_order_it == -orders[index])
        {
            this->_F_star[index] = std::conj(*F_minus_state_it);
            ++F_minus_order_it;
            ++F_minus_state_it;
        }
        if(Z_order_it != this->_orders.cend() && *Z_order_it == orders[index])
        {
            Z_new[index] = *Z_state_it;
            ++Z_order_it;
            ++Z_state_it;
        }
    }

    // … and finally update F̃^*(-0)
    this->_F_star[0] = std::conj(this->_F[0]);
    
    this->_Z = std::move(Z_new);
    this->_orders = std::move(orders);
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
    
    #pragma omp parallel for
    for(int order=0; order<this->_orders.size(); ++order)
    {
        this->_F[order] *= E.second;
        this->_F_star[order] *= E.second;
        this->_Z[order] *= E.first;
    }
    
    this->_Z[0] += 1.-E.first; // WARNING: assumes M0=1
}

void
Discrete
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
    
    #pragma omp parallel for
    for(int i=0; i<this->_orders.size(); ++i)
    {
        auto const k = this->_orders[i] * this->_bin_width;
        auto const D = operators::diffusion(
            this->species.get_D()[0].magnitude, duration.magnitude, 
            k.magnitude, delta_k);
        this->_F[i] *= std::get<0>(D);
        this->_F_star[i] *= std::get<1>(D);
        this->_Z[i] *= std::get<2>(D);
    }
}

void
Discrete
::off_resonance(Quantity const & duration)
{
    auto const angle = 
        duration * 2*M_PI*units::rad 
        * (this->delta_omega+this->species.get_delta_omega());
    if(angle.magnitude != 0)
    {
        auto const rotations = operators::phase_accumulation(angle.magnitude);
        
        #pragma omp parallel for
        for(int order=0; order<this->_orders.size(); ++order)
        {
            this->_states[0+3*order] *= rotations.first;
            this->_states[1+3*order] *= rotations.second;
            // Z̃ states are unaffected
        }
    }
}

Quantity const & 
Discrete
::bin_width() const
{
    return this->_bin_width;
}

}

}
