#include "Discrete.h"

#include <algorithm>
#include <set>
#include <vector>

#include "sycomore/epg/operators.h"
#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

Discrete
::Discrete(
    Species const & species, Magnetization const & initial_magnetization, 
    Quantity bin_width)
: species(species), _bin_width(bin_width), _states(3, 0)
{
    // Store magnetization as lines of F̃, F̃^*_, Z̃
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_states[0] = std::sqrt(2)*magnetization.p;
    this->_states[1] = std::sqrt(2)*magnetization.m;
    this->_states[2] = magnetization.z;
    
    this->_orders.push_back(0);
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
    return {
        this->_states[3*order],
        this->_states[3*order+1],
        this->_states[3*order+2]};
}

std::vector<Complex> const &
Discrete
::states() const
{
    return this->_states;
}

Complex const &
Discrete
::echo() const
{
    return this->_states[0];
}

void
Discrete
::apply_pulse(Quantity angle, Quantity phase)
{
    auto const T = operators::pulse(angle, phase);
    
    #pragma omp parallel for
    for(int order=0; order<this->_orders.size(); ++order)
    {
        std::vector<Complex> result(3, 0);
        for(int r=0; r<3; ++r)
        {
            for(int c=0; c<3; ++c)
            {
                result[r] += T[3*r+c] * this->_states[3*order+c];
            }
        }
        
        std::memcpy(
            this->_states.data()+3*order, result.data(), 3*sizeof(Complex));
    }
}

void
Discrete
::apply_time_interval(
    Quantity const & duration, Quantity const & gradient, Real threshold)
{
    if(duration.magnitude == 0)
    {
        return;
    }
    
    this->apply_relaxation(duration);
    this->apply_diffusion(duration, gradient);
    this->apply_gradient(duration, gradient);
    
    if(threshold > 0)
    {
        // TODO
        // abs = numpy.absolute(self.magnetization["v"])
        // keep = numpy.any(abs>threshold, axis=1)
        // keep[0] = True
        // self.magnetization = self.magnetization[keep]
    }
}

void
Discrete
::apply_gradient(Quantity const & duration, Quantity const & gradient)
{
    // This assumes a constant gradient in the integral: 
    // k(t) = γ ∫_0^t G(t') dt' = γ⋅t⋅G
    auto const delta_k = int64_t(std::round(
        (sycomore::gamma*gradient*duration / this->_bin_width).magnitude));
    
    if(delta_k == 0)
    {
        return;
    }
    
    // Unfold the F̃-states
    std::vector<int64_t> F_orders(2*this->_orders.size()-1, 0);
    std::vector<Complex> F(2*this->_orders.size()-1, 0);
    // F̃^+ on the right side
    for(int64_t i=0; i<this->_orders.size(); ++i)
    {
        F_orders[this->_orders.size()-1 + i] = this->_orders[i];
        F[this->_orders.size()-1 + i] = this->_states[3*i];
    }
    // F̃^{-*} on the left side, reversed order
    for(int64_t i=0; i<this->_orders.size(); ++i)
    {
        F_orders[this->_orders.size()-1 - i] = -this->_orders[i];
        F[this->_orders.size()-1 - i] = std::conj(this->_states[3*i+1]);
    }
    
    // Shift according to Δk
    for(auto && k: F_orders)
    {
        k += delta_k;
    }
    
    // Fold the F̃-states
    auto const Z_orders = this->_orders;
    std::set<int64_t> orders;
    for(auto && order: F_orders)
    {
        orders.insert(std::abs(order));
    }
    orders.insert(Z_orders.begin(), Z_orders.end());
    
    std::vector<Complex> magnetization(3*orders.size(), 0);
    
    uint destination_index=0;
    for(auto && order: orders)
    {
        auto F_orders_it = std::lower_bound(
            F_orders.begin(), F_orders.end(), order);
        if(F_orders_it != F_orders.end() && *F_orders_it <= order)
        {
            auto source_index = F_orders_it-F_orders.begin();
            magnetization[3*destination_index+0] = F[source_index];
        }
        
        F_orders_it = std::lower_bound(
            F_orders.begin(), F_orders.end(), -order);
        if(F_orders_it != F_orders.end() && *F_orders_it <= -order)
        {
            auto source_index = F_orders_it-F_orders.begin();
            magnetization[3*destination_index+1] = std::conj(F[source_index]);
        }
        
        auto const Z_orders_it = std::lower_bound(
            Z_orders.begin(), Z_orders.end(), order);
        if(Z_orders_it != Z_orders.end() && *Z_orders_it == order)
        {
            auto const source_index = Z_orders_it-Z_orders.begin();
            magnetization[3*destination_index+2] = this->_states[3*source_index+2];
        }
        
        ++destination_index;
    }
    magnetization[1] = std::conj(magnetization[0]);
    
    this->_states = magnetization;
    this->_orders = {orders.begin(), orders.end()};
}

void
Discrete
::apply_relaxation(Quantity const & duration)
{
    if(this->species.get_R1().magnitude == 0 && this->species.get_R2().magnitude == 0)
    {
        return;
    }
    
    auto const E = operators::relaxation(this->species, duration);
    
    #pragma omp parallel for
    for(int order=0; order<this->_orders.size(); ++order)
    {
        this->_states[0+3*order] *= E.second;
        this->_states[1+3*order] *= E.second;
        this->_states[2+3*order] *= E.first;
    }
    
    this->_states[2] += 1.-E.first; // WARNING: assumes M0=1
}

void
Discrete
::apply_diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(this->species.get_D().magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = sycomore::gamma*gradient*duration;
    
    #pragma omp parallel for
    for(int i=0; i<this->_orders.size(); ++i)
    {
        auto const k = this->_orders[i] * this->_bin_width;
        auto const D = operators::diffusion(this->species, duration, k, delta_k);
        this->_states[0+3*i] *= std::get<0>(D);
        this->_states[1+3*i] *= std::get<1>(D);
        this->_states[2+3*i] *= std::get<2>(D);
    }
}

}

}
