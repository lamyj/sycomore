#include "Discrete3D.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

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
: species(species), _bin_width(bin_width), _states()
{
    auto const M = as_complex_magnetization(initial_magnetization);
    this->_states[{0,0,0}] = {std::sqrt(2)*M.p, std::sqrt(2)*M.m, M.z};
}

std::vector<Complex>
Discrete3D
::state(Array<Quantity> const & order) const
{
    Array<int64_t> bin(order.size());
    std::transform(
        order.begin(), order.end(), bin.begin(),
        [&](Quantity const & q) { return std::round((q/this->_bin_width).magnitude); });
    return this->_states.at(bin);
}

std::unordered_map<Array<Quantity>, std::vector<Complex>>
Discrete3D
::states() const
{
    std::unordered_map<Array<Quantity>, std::vector<Complex>> states;
    for(auto && state: this->_states)
    {
        states[state.first*this->_bin_width] = state.second;
    }
    return states;
}

Complex const &
Discrete3D
::echo() const
{
    return this->_states.at({0,0,0})[0];
}

void
Discrete3D
::apply_pulse(Quantity angle, Quantity phase)
{
    auto const T = operators::pulse(angle, phase);

    // FIXME: OpenMP
    // Not detected as a loop construct
//    for(auto && state: this->_states)
    // Condition of OpenMP for loop must be a relational comparison ('<', '<=', '>', or '>=')
    // unordered_map does not have comparable iterators
//    for(auto it=this->_states.begin(); it!=this->_states.end(); ++it)

#if defined(_OPENMP)
    // Dirty version: build a vector of iterators and iterate on it.
    std::vector<decltype(this->_states.begin())> iterators;
    iterators.reserve(this->_states.size());
    for(auto it=this->_states.begin(); it!=this->_states.end(); ++it)
    {
        iterators.push_back(it);
    }
    #pragma omp parallel for
    for(auto it=iterators.begin(); it<iterators.end(); ++it)
    {
        auto && state = **it;
#else
    for(auto && state: this->_states)
    {
#endif
        std::vector<Complex> result(3, 0);
        for(int r=0; r<3; ++r)
        {
            for(int c=0; c<3; ++c)
            {
                result[r] += T[3*r+c] * state.second[c];
            }
        }
        state.second = std::move(result);
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
        decltype(this->_states) states;
        for(auto && state: this->_states)
        {
            auto && population = state.second;
            auto const magnitude_squared =
                std::pow(std::abs(population[0]), 2)
                +std::pow(std::abs(population[1]), 2)
                +std::pow(std::abs(population[2]), 2);
            if(magnitude_squared >= threshold_squared)
            {
                states.insert(state);
            }
        }

        this->_states = std::move(states);
    }
}

void
Discrete3D
::shift(Quantity const & duration, Array<Quantity> const & gradient)
{
    // This assumes a constant gradient in the integral:
    // k(t) = γ ∫_0^t G(t') dt' = γ⋅t⋅G
    auto const delta_k_q = sycomore::gamma*gradient*duration;
    Array<int64_t> delta_k(delta_k_q.size());
    std::transform(
        delta_k_q.begin(), delta_k_q.end(), delta_k.begin(),
        [&](Quantity const & q) {
            return (int64_t)(std::round((q/this->_bin_width).magnitude)); });

    decltype(this->_states) states(this->_states.bucket_count());
    states[{0,0,0}] = {0,0,0};
    for(auto state: this->_states)
    {
        auto && order = state.first;
        auto && population = state.second;

        // Z̃ orders are not shifted
        {
            auto && it = states.insert({order, {0,0,0}}).first;
            it->second[2] = population[2];
        }

        // Shift F̃_{+}(k)
        {
            auto const new_order = order+delta_k;
            if(new_order[0] >= 0)
            {
                // The state stays in the same half-space after shifting.
                // Store as F̃_{+}(k+delta_k)
                auto && it = states.insert({new_order, {0,0,0}}).first;
                it->second[0] = population[0];
            }
            else
            {
                // The state moves to the other half-space after shifting.
                // Store as F̃^{*}_{-}(-k+delta_k)
                auto && it = states.insert({-new_order, {0,0,0}}).first;
                it->second[1] = std::conj(population[0]);
            }
        }

        // Shift F̃^{*}_{-}(-k)
        {
            auto const new_order = -order+delta_k;
            if(new_order[0] < 0)
            {
                // The state stays in the same half-space after shifting.
                // Store as F̃^{*}_{-}(-k+delta_k)
                auto && it = states.insert({-new_order, {0,0,0}}).first;
                it->second[1] = population[1];
            }
            else
            {
                // The state moves to the other half-space after shifting.
                // Store as F̃^_{+}(k-delta_k)
                auto && it = states.insert({new_order, {0,0,0}}).first;
                it->second[0] = std::conj(population[1]);
            }
        }

    }

    this->_states = std::move(states);
}

void
Discrete3D
::relaxation(Quantity const & duration)
{
    auto const E = operators::relaxation(this->species, duration);
    // FIXME: OpenMP
    for(auto && state: this->_states)
    {
        state.second[0] *= E.second;
        state.second[1] *= E.second;
        state.second[2] *= E.first;
    }
    this->_states[{0,0,0}][2] += 1.-E.first; // WARNING: assumes M0=1
}

void
Discrete3D
::diffusion(Quantity const & duration, Array<Quantity> const & gradient)
{
    auto const delta_k = sycomore::gamma*gradient*duration;

    // FIXME: OpenMP
    for(auto && state: this->_states)
    {
        auto const k1 = state.first * this->_bin_width;
        auto const k2_plus = k1+delta_k;
        auto const k2_minus = -k1+delta_k;

        // Weigel 2010, eq. 26
        Quantity b_L[3][3];
        for(std::size_t m=0; m<3; ++m)
        {
            for(std::size_t n=0; n<3; ++n)
            {
                b_L[m][n] = k1[m]*k1[n]*duration;
            }
        }

        // Weigel 2010, eq. 30
        Quantity b_T_plus[3][3];
        Quantity b_T_minus[3][3];
        for(std::size_t m=0; m<3; ++m)
        {
            for(std::size_t n=0; n<3; ++n)
            {
                b_T_plus[m][n] =
                    b_L[m][n]
                    + 0.5*k1[m]*(k2_plus[n]-k1[n])*duration
                    + 0.5*k1[n]*(k2_plus[m]-k1[m])*duration
                    + 1./3.*(k2_plus[m]-k1[m])*(k2_plus[n]-k1[n])*duration;
                b_T_minus[m][n] =
                    b_L[m][n]
                    + 0.5*k1[m]*(k2_minus[n]-k1[n])*duration
                    + 0.5*k1[n]*(k2_minus[m]-k1[m])*duration
                    + 1./3.*(k2_minus[m]-k1[m])*(k2_minus[n]-k1[n])*duration;
            }
        }

        // FIXME
        // this->species.D should be a tensor
        auto const D_T_plus = exp(
            (-(b_T_plus[0][0]+b_T_plus[1][1]+b_T_plus[2][2])*this->species.get_D()).magnitude);
        auto const D_T_minus = exp(
            (-(b_T_minus[0][0]+b_T_minus[1][1]+b_T_minus[2][2])*this->species.get_D()).magnitude);
        auto const D_L = exp(
            (-(b_L[0][0]+b_L[1][1]+b_L[2][2])*this->species.get_D()).magnitude);

        state.second[0] *= D_T_plus;
        state.second[1] *= D_T_plus;
        state.second[2] *= D_L;
    }
}

}

}
