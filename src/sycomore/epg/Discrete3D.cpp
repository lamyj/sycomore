#include "Discrete3D.h"

#include <algorithm>
#include <cmath>
#include <cstring>
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

    // Unfold the F̃ states
    std::vector<decltype(this->_orders)::value_type> F_orders((2*this->size()-1)*3);
    std::vector<decltype(this->_states)::value_type> F_states(2*this->size()-1);
    for(std::size_t i=0; i<this->size(); ++i)
    {
        // F̃^+ on the right side.
        F_orders[3*(this->size()-1 + i)+0] = this->_orders[3*i+0];
        F_orders[3*(this->size()-1 + i)+1] = this->_orders[3*i+1];
        F_orders[3*(this->size()-1 + i)+2] = this->_orders[3*i+2];

        F_states[this->size()-1 + i] = this->_states[3*i+0];

        // F̃^{-*} on the left side, reversed order.
        F_orders[3*(this->size()-1 - i)+0] = -this->_orders[3*i+0];
        F_orders[3*(this->size()-1 - i)+1] = -this->_orders[3*i+1];
        F_orders[3*(this->size()-1 - i)+2] = -this->_orders[3*i+2];

        F_states[this->size()-1 - i] = std::conj(this->_states[3*i+1]);
    }

    // Shift according to Δk
    for(std::size_t i=0; i<F_states.size(); ++i)
    {
        for(std::size_t j=0; j<3; ++j)
        {
            F_orders[3*i+j] += delta_k[j];
        }
    }

    // Find the location of the negative and positive orders.
    // NOTE The non-const iterator is required in order to build F_plus_order
    // as a view in the loop
    std::pair<decltype(F_orders.begin()), decltype(F_states.begin())> F_plus_it(
        F_orders.end(), F_states.end());
    for(auto it = F_orders.begin(); it != F_orders.end(); it += 3)
    {
        if(*it >= 0)
        {
            F_plus_it.first = it;
            F_plus_it.second = F_states.begin()+(it-F_orders.begin())/3;
            break;
        }
    }
    // NOTE Since we use -order, we cannot build a view
    std::pair<decltype(F_orders.begin()), decltype(F_states.begin())> F_minus_it(
        F_orders.begin()-3, F_states.begin()-1);
    for(auto it = F_orders.end()-3; it >= F_orders.begin(); it -= 3)
    {
        if(*it < 0)
        {
            F_minus_it.first = it;
            F_minus_it.second = F_states.begin()+(it-F_orders.begin())/3;
            break;
        }
    }

    // NOTE The non-const iterator is required in order to build Z_order as a
    // view in the loop
    std::pair<decltype(this->_orders.begin()), decltype(this->_states.begin())> Z_it(
        this->_orders.begin(), this->_states.begin());

    decltype(this->_orders) orders; orders.reserve(3*(this->size()+F_states.size()/2));
    decltype(this->_states) states; states.reserve(3*(this->size()+F_states.size()/2));
    while(!(
        F_plus_it.first == F_orders.end()
        && F_minus_it.first < F_orders.begin()
        && Z_it.first == this->_orders.end()))
    {
        // Create the orders from the iterators
        auto const F_plus_order =
            F_plus_it.first != F_orders.end()
            ? Bin(&*F_plus_it.first, 3) : Bin();
        auto const F_minus_order =
            F_minus_it.first >= F_orders.begin()
            ? -Bin{*(F_minus_it.first+0), *(F_minus_it.first+1), *(F_minus_it.first+2)}
            : Bin();
        auto const Z_order =
            Z_it.first != this->_orders.end() ? Bin(&*Z_it.first, 3) : Bin();

        // Find the orders which have the minimum value.
        auto const do_F_plus =
            !F_plus_order.empty()
            && (F_minus_order.empty() || F_plus_order <= F_minus_order)
            && (Z_order.empty() || F_plus_order <= Z_order);
        auto const do_F_minus =
            !F_minus_order.empty()
            && (F_plus_order.empty() || F_minus_order <= F_plus_order)
            && (Z_order.empty() || F_minus_order <= Z_order);
        auto const do_Z =
            !Z_order.empty()
            && (F_plus_order.empty() || Z_order <= F_plus_order)
            && (F_minus_order.empty() || Z_order <= F_minus_order);

        auto && order =
            do_F_plus ? F_plus_order : (do_F_minus ? F_minus_order: Z_order);

        for(auto && x: order)
        {
            orders.push_back(x);
        }

        if(do_F_plus)
        {
            states.push_back(*F_plus_it.second);
            F_plus_it.first += 3;
            ++F_plus_it.second;
        }
        else
        {
            states.push_back(0);
        }

        if(do_F_minus)
        {
            states.push_back(std::conj(*F_minus_it.second));
            F_minus_it.first -= 3;
            --F_minus_it.second;
        }
        else
        {
            states.push_back(0);
        }

        if(do_Z)
        {
            states.push_back(*(Z_it.second+2));
            Z_it.first += 3;
            Z_it.second += 3;
        }
        else
        {
            states.push_back(0);
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
    for(std::size_t order=0; order<this->size(); ++order)
    {
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

}

}
