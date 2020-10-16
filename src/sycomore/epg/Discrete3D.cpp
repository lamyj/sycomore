#include "Discrete3D.h"

#include <algorithm>
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
    this->_zero_it = 0;
}

std::size_t
Discrete3D
::size() const
{
    return this->_orders.size() / 3;
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

    using namespace sycomore::units;

    Bin bin(order.size());
    for(std::size_t i=0; i<bin.size(); ++i)
    {
        bin[i] = order[i]/this->_bin_width;
    }
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
    return this->_F[this->_zero_it];
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
        decltype(this->_orders) orders; orders.reserve(this->_orders.size());
        decltype(this->_F) F; F.reserve(this->_F.size());
        decltype(this->_F_star) F_star; F_star.reserve(this->_F_star.size());
        decltype(this->_Z) Z; Z.reserve(this->_Z.size());

        for(std::size_t index=0; index<this->size(); ++index)
        {
            Bin const order(this->_orders.data()+3*index, 3);
            
            auto const magnitude_squared =
                std::pow(std::abs(this->_F[index]), 2)
                +std::pow(std::abs(this->_F_star[index]), 2)
                +std::pow(std::abs(this->_Z[index]), 2);
            if(magnitude_squared >= threshold_squared)
            {
                std::copy(order.begin(), order.end(), std::back_inserter(orders));
                F.push_back(this->_F[index]);
                F_star.push_back(this->_F_star[index]);
                Z.push_back(this->_Z[index]);
            }
        }

        this->_orders = std::move(orders);
        this->_F = std::move(F);
        this->_F_star = std::move(F_star);
        this->_Z = std::move(Z);
        
        // Update the iterator pointing to the echo magnetization.
        for(std::size_t i=0; i<this->size(); ++i)
        {
            Bin const order(this->_orders.data()+3*i, 3);
            if(order == Array<int64_t>{0,0,0})
            {
                this->_zero_it = i;
                break;
            }
        }
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
    // WARNING: this does not work if delta_k is zero
    Bin delta_k(gradient.size());
    bool all_zero = true;
    for(std::size_t i=0; i<delta_k.size(); ++i)
    {
        delta_k[i] = std::round(
            sycomore::gamma.magnitude*gradient[i].magnitude*duration.magnitude
            / this->_bin_width.magnitude);
        if(delta_k[i] != 0)
        {
            all_zero = false;
        }
    }
    if(all_zero)
    {
        return;
    }
    
    // std::cout << "Before shift\n";
    // for(int i=0; i<this->size(); ++i)
    // {
    //     std::cout 
    //         << this->_orders[3*i] << " " << this->_orders[3*i+1] << " " << this->_orders[3*i+2] << ": " 
    //         << this->_F[i] << " " << this->_F_star[i] << " " << this->_Z[i] << "\n"; 
    // }
    // std::cout << "\n";
    
    /*
     * The shift operator is the only operator which will create new orders, and
     * which will move population from one state to another. Other operators 
     * only change the population of each state independently. The shift 
     * operator must then maintain the relative order of the states. Any order
     * can be used, as long as F states (stored as a conjugate pair) can be
     * unfolded and re-folded while keeping the order. In this implementation,
     * we use the lexicographic order (F_x < F_y < F_z), and the orders vector
     * will only contain orders in the upper quadrant (F_x <= F_y <= F_z).
     * 
     * The shift operator starts by unfolding the conjugate pairs of F, the 
     * same operation is applied to the vector of orders. These two unfolded 
     * vectors needs not be sorted.
     *
     * All orders are shifted by delta_k. For an arbitrary delta_k, this will
     * break the order of sorted vector, hence the non-sorted requirement above.
     *
     * The unfolded F-orders/F-states and the non-shifted Z-order/Z-states are
     * then inserted in new orders/states vectors.
     */

    // Unfold the F̃ states
    decltype(this->_orders) F_orders; F_orders.reserve((2*this->size())*3);
    decltype(this->_F) F_states; F_states.reserve(2*this->size());
    for(std::size_t i=0; i<this->size(); ++i)
    {
        if(this->_F[i] != 0.)
        {
            F_orders.push_back(this->_orders[3*i+0]);
            F_orders.push_back(this->_orders[3*i+1]);
            F_orders.push_back(this->_orders[3*i+2]);
            F_states.push_back(this->_F[i]);
        }
        
        if(i != this->_zero_it && this->_F_star[i] != 0.)
        {
            F_orders.push_back(-this->_orders[3*i+0]);
            F_orders.push_back(-this->_orders[3*i+1]);
            F_orders.push_back(-this->_orders[3*i+2]);
            F_states.push_back(std::conj(this->_F_star[i]));
        }
    }
    
    // std::cout << "Unfolded\n";
    // for(int i=0; i<F_states.size(); ++i)
    // {
    //     std::cout 
    //         << F_orders[3*i] << " " << F_orders[3*i+1] << " " << F_orders[3*i+2] << ": " 
    //         << F_states[i] << "\n"; 
    // }
    // std::cout << "\n";
    
    // Shift according to Δk
    for(std::size_t i=0; i<F_states.size(); ++i)
    {
        for(std::size_t j=0; j<3; ++j)
        {
            F_orders[3*i+j] += delta_k[j];
        }
    }
    
    // New orders, pre-filled with order 0. Cannot be modified in-place,
    // we need the old orders when inserting the Z orders/states.
    decltype(this->_orders) orders; 
    orders.reserve(this->_orders.size());
    orders.push_back(0); orders.push_back(0); orders.push_back(0);
    // New F states, can be modified in-place.
    this->_F.clear();
    this->_F.push_back(0);
    // Same for F* states.
    this->_F_star.clear();
    this->_F_star.push_back(0);
    // New Z states. Cannot be modified in-place, cf. remark on old orders. 
    decltype(this->_Z) Z; 
    Z.reserve(this->_Z.size());
    Z.push_back(0);
    // Mapping between a normalized (i.e. folded) order and its location in the
    // states vectors.
    std::map<Bin, std::size_t> locations;
    locations[{0,0,0}] = 0;
    for(std::size_t i=0; i != F_states.size(); ++i)
    {
        auto & value = F_states[i];
        
        Bin order{F_orders[3*i+0], F_orders[3*i+1], F_orders[3*i+2]};
        auto array = &this->_F;
        if(order[0] < 0)
        {
            order *= -1;
            value = std::conj(value);
            array = &this->_F_star;
        }
        
        auto locations_it = locations.find(order);
        if(locations_it == locations.end())
        {
            this->_F.push_back(0);
            this->_F_star.push_back(0);
            Z.push_back(0);
            std::copy(order.begin(), order.end(), std::back_inserter(orders));
            locations_it = 
                locations.emplace(std::move(order), this->_F.size()-1).first;
        }
        (*array)[locations_it->second] = value;
    }
    for(std::size_t i=0; i != this->size(); ++i)
    {
        auto const value = this->_Z[i];
        if(value == 0.)
        {
            continue;
        }
        
        Bin const order{this->_orders.data()+3*i, 3};
        auto locations_it = locations.find(order);
        if(locations_it == locations.end())
        {
            this->_F.push_back(0);
            this->_F_star.push_back(0);
            Z.push_back(0);
            std::copy(order.begin(), order.end(), std::back_inserter(orders));
            locations_it = locations.emplace(order, this->_F.size()-1).first;
        }
        Z[locations_it->second] = value;
    }
    
    this->_orders = std::move(orders);
    this->_Z = std::move(Z);

    // std::cout << "Folded Model\n";
    // for(int i=0; i<this->size(); ++i)
    // {
    //     std::cout 
    //         << this->_orders[3*i] << " " << this->_orders[3*i+1] << " " << this->_orders[3*i+2] << ": " 
    //         << this->_F[i] << " " << this->_F_star[i] << " " << this->_Z[i] << "\n"; 
    // }
    // std::cout << "\n";

    // Make sure the iterator pointing to the echo magnetization is at the 
    // correct position, update the conjugate magnetization of the echo 
    // magnetization.
    this->_zero_it = 0;
    this->_F_star[this->_zero_it] = std::conj(this->_F[this->_zero_it]);
    
    // std::cout << "After shift\n";
    // for(int i=0; i<this->size(); ++i)
    // {
    //     std::cout 
    //         << this->_orders[3*i] << " " << this->_orders[3*i+1] << " " << this->_orders[3*i+2] << ": " 
    //         << this->_F[i] << " " << this->_F_star[i] << " " << this->_Z[i] << "\n"; 
    // }
    // std::cout << "k=0: " << this->_zero_it << "\n";
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
    
    this->_Z[this->_zero_it] += 1.-E.first; // WARNING: assumes M0=1
}

void
Discrete3D
::diffusion(Quantity const & duration, Array<Quantity> const & gradient)
{
    using namespace sycomore::units;

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
    
    auto const tau = duration.convert_to(s);
    auto const bin_width = this->_bin_width.convert_to(rad/m);
    
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
    
    std::vector<Real> const delta_k{
        (sycomore::gamma*gradient[0]*duration).convert_to(rad/m),
        (sycomore::gamma*gradient[1]*duration).convert_to(rad/m),
        (sycomore::gamma*gradient[2]*duration).convert_to(rad/m)
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
                this->species.get_D()[3*m+n].convert_to(units::m*units::m/units::s), 
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
        duration * 2*M_PI*units::rad 
        * (this->delta_omega+this->species.get_delta_omega());
    if(angle.magnitude != 0)
    {
        auto const rotations = operators::phase_accumulation(angle.magnitude);
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
