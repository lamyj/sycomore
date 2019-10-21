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
#include "sycomore/sycomore.h"

namespace sycomore
{

namespace epg
{

Discrete3D
::Discrete3D(
    Species const & species, Magnetization const & initial_magnetization,
    Quantity bin_width)
: species(species), _bin_width(bin_width)
{
    auto const M = as_complex_magnetization(initial_magnetization);
    this->_orders = {0,0,0};
    this->_states = {std::sqrt(2)*M.p, std::sqrt(2)*M.m, M.z};
    this->_zero_it = this->_states.begin();
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
::state(Order const & order) const
{
    if(order.size() != 3)
    {
        std::ostringstream message;
        message<< "Order must have 3 elements, not " << order.size();
        throw std::runtime_error(message.str());
    }

    using namespace sycomore::units;

    auto const bin_width = this->_bin_width.convert_to(rad/m);
    Bin bin(order.size());
    for(std::size_t i=0; i<bin.size(); ++i)
    {
        bin[i] = order[i].convert_to(rad/m)/bin_width;
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

    auto const index = it-this->_orders.begin();
    return {
        this->_states[index+0], this->_states[index+1], this->_states[index+2]};
}

std::vector<Complex> const &
Discrete3D
::states() const
{
    return this->_states;
}

Complex const &
Discrete3D
::echo() const
{
    return *this->_zero_it;
}

void
Discrete3D
::apply_pulse(Quantity angle, Quantity phase)
{
    auto const T = operators::pulse(angle, phase);

    #pragma omp parallel for
    for(int order=0; order<this->size(); ++order)
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
Discrete3D
::apply_time_interval(
    Quantity const & duration, Array<Quantity> const & gradient, Real threshold)
{
    if(duration.magnitude == 0)
    {
        return;
    }

    this->relaxation(duration);
    this->diffusion(duration, gradient);
    this->shift(duration, gradient);

    if(threshold > 0)
    {
        auto const threshold_squared = std::pow(threshold, 2);
        decltype(this->_orders) orders; orders.reserve(this->_orders.size());
        decltype(this->_states) states; states.reserve(this->_states.size());

        for(std::size_t i=0; i<this->size(); ++i)
        {
            Bin const order(this->_orders.data()+3*i, 3);
            Array<Complex> const state(this->_states.data()+3*i, 3);
            auto const magnitude_squared =
                std::pow(std::abs(state[0]), 2)
                +std::pow(std::abs(state[1]), 2)
                +std::pow(std::abs(state[2]), 2);
            if(magnitude_squared >= threshold_squared || order == Bin{0,0,0})
            {
                std::copy(order.begin(), order.end(), std::back_inserter(orders));
                std::copy(state.begin(), state.end(), std::back_inserter(states));
            }
        }

        this->_orders = std::move(orders);
        this->_states = std::move(states);

        // Update the iterator pointing to the echo magnetization.
        for(std::size_t i=0; i<this->size(); ++i)
        {
            Bin const order(this->_orders.data()+3*i, 3);
            if(order == Array<int64_t>{0,0,0})
            {
                this->_zero_it = this->_states.begin()+3*i;
                break;
            }
        }
    }
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
        delta_k[i] = std::round((delta_k_q[i]/this->_bin_width).magnitude);
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
    //         << this->_states[3*i] << " " << this->_states[3*i+1] << " " << this->_states[3*i+2] << "\n"; 
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
     * The shift operator starts by unfolding the conjugate pairs of F states,
     * creating a vector of size 2*N-1 from a vector of size N (the order 0 must
     * not be unfolded). The same operation is applied to the vector of orders.
     * These two unfolded vectors needs not be sorted.
     *
     * All orders are shifted by delta_k. For an arbitrary delta_k, this will
     * break the order of sorted vector, hence the non-sorted requirement above.
     *
     * A map of the folded orders is then built, and split between this->_orders
     * and this->_states.
     */

    // Unfold the F̃ states
    std::vector<decltype(this->_orders)::value_type> F_orders((2*this->size())*3);
    std::vector<decltype(this->_states)::value_type> F_states(2*this->size());
    for(std::size_t i=0; i<this->size(); ++i)
    {
        // F̃^+ on the left side.
        F_orders[3*i+0] = this->_orders[3*i+0];
        F_orders[3*i+1] = this->_orders[3*i+1];
        F_orders[3*i+2] = this->_orders[3*i+2];
    
        F_states[i] = this->_states[3*i+0];
    
        // F̃^{-*} on the right side.
        F_orders[3*(this->size()+i)+0] = -this->_orders[3*i+0];
        F_orders[3*(this->size()+i)+1] = -this->_orders[3*i+1];
        F_orders[3*(this->size()+i)+2] = -this->_orders[3*i+2];
        F_states[this->size()+i] = std::conj(this->_states[3*i+1]);
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
        Bin order;
        unsigned int location;
        decltype(this->_states)::value_type value;
        if(*F_orders_it <= *(F_orders_it+1) && *(F_orders_it+1) <= *(F_orders_it+2))
        {
            // Upper quadrant: use as is
            order = Bin(&*F_orders_it, 3);
            location = 0;
            value = *F_states_it;
        }
        else
        {
            // Lower quadrant: store conjugate
            order = -Bin{*(F_orders_it+0), *(F_orders_it+1), *(F_orders_it+2)};
            location = 1;
            value = std::conj(*F_states_it);
        }
    
        // Create empty value if needed
        auto it = folded.insert({order, {0,0,0}}).first;
    
        it->second[location] = value;
    
        F_orders_it += 3;
        ++F_states_it;
    }
    auto && Z_orders_it = this->_orders.begin();
    auto && Z_states_it = this->_states.begin();
    while(Z_orders_it != this->_orders.end())
    {
        Bin const order(&*Z_orders_it, 3);
        unsigned int const location=2;
        decltype(this->_states)::value_type const value=*(Z_states_it+2);
    
        // Create empty value if needed
        auto it = folded.insert({order, {0,0,0}}).first;
    
        it->second[location] = value;
    
        Z_orders_it += 3;
        Z_states_it += 3;
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
    this->_states.resize(3*folded.size());
    auto orders_it = this->_orders.begin();
    auto states_it = this->_states.begin();
    for(auto && item: folded)
    {
        orders_it = std::copy(item.first.begin(), item.first.end(), orders_it);
        states_it = std::copy(item.second.begin(), item.second.end(), states_it);
    }
    
    // std::cout << "Folded Model\n";
    // for(int i=0; i<this->size(); ++i)
    // {
    //     std::cout 
    //         << this->_orders[3*i] << " " << this->_orders[3*i+1] << " " << this->_orders[3*i+2] << ": " 
    //         << this->_states[3*i] << " " << this->_states[3*i+1] << " " << this->_states[3*i+2] << "\n"; 
    // }
    // std::cout << "\n";

    // Update the iterator pointing to the echo magnetization.
    for(std::size_t i=0; i<this->size(); ++i)
    {
        Bin const order(this->_orders.data()+3*i, 3);
        if(order == Array<int64_t>{0,0,0})
        {
            this->_zero_it = this->_states.begin()+3*i;
            break;
        }
    }
    // Update the conjugate magnetization of the echo magnetization
    this->_zero_it[1] = std::conj(this->_zero_it[0]);
}

void
Discrete3D
::relaxation(Quantity const & duration)
{
    if(this->species.get_R1().magnitude == 0 && this->species.get_R2().magnitude == 0)
    {
        return;
    }

    auto const E = operators::relaxation(this->species, duration);

    #pragma omp parallel for
    for(int order=0; order<this->size(); ++order)
    {
        this->_states[0+3*order] *= E.second;
        this->_states[1+3*order] *= E.second;
        this->_states[2+3*order] *= E.first;
    }

    this->_zero_it[2] += 1.-E.first; // WARNING: assumes M0=1
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
        D[i] = this->species.get_D()[i].convert_to(m*m/s);
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
                        + 1./3.*delta_k[m]*bin_width*delta_k[n]);
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

        this->_states[3*order+0] *= exp(-b_T_plus_D);
        this->_states[3*order+1] *= exp(-b_T_minus_D);
        this->_states[3*order+2] *= exp(-b_L_D);
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
