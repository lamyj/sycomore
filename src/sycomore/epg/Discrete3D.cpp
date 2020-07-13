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
    auto const T = operators::pulse(angle.magnitude, phase.magnitude);

    for(std::size_t order = 0; order < this->size(); ++order)
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
    auto const delta_k_q = sycomore::gamma*gradient*duration;
    Bin delta_k(delta_k_q.size());
    bool all_zero = true;
    for(std::size_t i=0; i<delta_k.size(); ++i)
    {
        delta_k[i] = std::round((delta_k_q[i]/this->_bin_width));
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
     * A map of the folded orders and the non-shifted Z-orders is then built.
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
    
    // Create the folded map for F and Z states, make sure the 0 state exists.
    std::map<Bin, State> folded({{{0,0,0}, {0,0,0}}});
    auto && F_orders_it = F_orders.begin();
    auto && F_states_it = F_states.begin();
    while(F_orders_it != F_orders.end())
    {
        Bin order{*(F_orders_it+0), *(F_orders_it+1), *(F_orders_it+2)};
        unsigned int location;
        decltype(F_states)::value_type value;
        if(order[0] >= 0)
        {
            // Right half-space: use as is
            location = 0;
            value = *F_states_it;
        }
        else
        {
            // Left half-space: use conjugate
            order *= -1;
            location = 1;
            value = std::conj(*F_states_it);
        }
    
        // Create empty value if needed
        auto it = folded.emplace(std::move(order), State{0,0,0}).first;
        it->second[location] = value;
    
        F_orders_it += 3;
        ++F_states_it;
    }
    auto && Z_orders_it = this->_orders.begin();
    auto && Z_states_it = this->_Z.begin();
    while(Z_orders_it != this->_orders.end())
    {
        auto const & value = *Z_states_it;
        if(value != 0.)
        {
            Bin const order(&*Z_orders_it, 3);
            unsigned int const location = 2;
            
            // Create empty value if needed
            auto it = folded.emplace(order, State{0,0,0}).first;
            it->second[location] = value;
        }
    
        Z_orders_it += 3;
        ++Z_states_it;
    }
    
    // std::cout << "Folded map\n";
    // for(auto & item: folded)
    // {
    //     std::cout 
    //         << item.first << ": "
    //         << item.second[0] << " " << item.second[1] << " " << item.second[2] << "\n";
    // }
    // std::cout << "\n";
    
    this->_orders.resize(3*folded.size());
    this->_F.resize(folded.size());
    this->_F_star.resize(folded.size());
    this->_Z.resize(folded.size());
    auto orders_it = this->_orders.begin();
    auto F_it = this->_F.begin();
    auto F_star_it = this->_F_star.begin();
    auto Z_it = this->_Z.begin();
    for(auto && item: folded)
    {
        orders_it = std::copy(item.first.begin(), item.first.end(), orders_it);
        *F_it = item.second[0];
        *F_star_it = item.second[1];
        *Z_it = item.second[2];
        
        ++F_it;
        ++F_star_it;
        ++Z_it;
    }
    
    // std::cout << "Folded Model\n";
    // for(int i=0; i<this->size(); ++i)
    // {
    //     std::cout 
    //         << this->_orders[3*i] << " " << this->_orders[3*i+1] << " " << this->_orders[3*i+2] << ": " 
    //         << this->_F[i] << " " << this->_F_star[i] << " " << this->_Z[i] << "\n"; 
    // }
    // std::cout << "\n";

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
    // Update the conjugate magnetization of the echo magnetization
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

    for(int order=0; order<this->size(); ++order)
    {
        this->_F[order] *= E.second;
        this->_F_star[order] *= E.second;
        this->_Z[order] *= E.first;
    }

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

    std::vector<Real> const delta_k{
        (sycomore::gamma*gradient[0]*duration).convert_to(rad/m),
        (sycomore::gamma*gradient[1]*duration).convert_to(rad/m),
        (sycomore::gamma*gradient[2]*duration).convert_to(rad/m)
    };

    Real D[9];
    for(std::size_t i=0; i<9; ++i)
    {
        D[i] = this->species.get_D()[i].magnitude;
    }
    auto const tau = duration.convert_to(s);
    auto const bin_width = this->_bin_width.convert_to(rad/m);

    #pragma omp parallel for
#ifdef _WIN32
    // WARNING: only signed integer types in OpenMP loops on Windows
    for(int order=0; order<this->size(); ++order)
    {
#else
    for(std::size_t order=0; order<this->size(); ++order)
    {
#endif
        // NOTE: this is not really k1, but rather k1/bin_width. However,
        // not declaring a new variable to store k1 is much faster.
        auto const k1=this->_orders.data()+3*order;

        // Weigel 2010, eq. 26
        Real b_L[9];
        int i=0;
        for(std::size_t m=0; m<3; ++m)
        {
            for(std::size_t n=0; n<3; ++n, ++i)
            {
                b_L[i] = (k1[m]*bin_width)*(k1[n]*bin_width)*tau;
            }
        }

        // Weigel 2010, eq. 30. Use delta_k instead of k2-k1, and factorize tau.
        Real b_T_plus[9];
        Real b_T_minus[9];
        i=0;
        for(std::size_t m=0; m<3; ++m)
        {
            for(std::size_t n=0; n<3; ++n, ++i)
            {
                b_T_plus[i] =
                    b_L[i] + tau * (
                          1./2.*(k1[m]*bin_width)*delta_k[n]
                        + 1./2.*(k1[n]*bin_width)*delta_k[m]
                        + 1./3.*delta_k[m]*delta_k[n]);
                b_T_minus[i] =
                    b_L[i] + tau * (
                          1./2.*(-k1[m]*bin_width)*delta_k[n]
                        + 1./2.*(-k1[n]*bin_width)*delta_k[m]
                        + 1./3.*delta_k[m]*delta_k[n]);
            }
        }

        // NOTE eq. 32 and 33 in Weigel 2010 use the trace of a product of two
        // matrices. This is equal to the sum of entry-wise product of elements.
        Real b_T_plus_D=0, b_T_minus_D=0, b_L_D=0;
        for(std::size_t i=0; i<9; ++i)
        {
            b_T_plus_D += b_T_plus[i]*D[i];
            b_T_minus_D += b_T_minus[i]*D[i];
            b_L_D += b_L[i]*D[i];
        }

        this->_F[order] *= exp(-b_T_plus_D);
        this->_F_star[order] *= exp(-b_T_minus_D);
        this->_Z[order] *= exp(-b_L_D);
    }
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
        
        for(int order=0; order<this->size(); ++order)
        {
            this->_F[order] *= rotations.first;
            this->_F_star[order] *= rotations.second;
            // Z̃ states are unaffected
        }
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
