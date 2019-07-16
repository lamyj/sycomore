#include "Regular.h"

#include <cstring>
#include <vector>

#include "sycomore/Array.h"
#include "sycomore/epg/operators.h"
#include "sycomore/Grid.h"
#include "sycomore/GridScanner.h"
#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

Regular
::Regular(
    Species const & species, Magnetization const & initial_magnetization, 
    unsigned int initial_size)
: species(species), _states(3*initial_size, 0)
{
    // Store magnetization as lines of F̃_k, F̃^*_{-k}, Z̃_k
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_states[0] = std::sqrt(2)*magnetization.p;
    this->_states[1] = std::sqrt(2)*magnetization.m;
    this->_states[2] = magnetization.z;
    
    this->_states_count = 1;
}

std::size_t const 
Regular
::states_count() const
{
    return this->_states_count;
}

std::vector<Complex>
Regular
::state(std::size_t order) const
{
    return {
        this->_states[3*order], 
        this->_states[3*order+1], 
        this->_states[3*order+2]};
}

std::vector<Complex>
Regular
::states() const
{
    return std::vector<Complex>{
        this->_states.begin(), this->_states.begin()+3*this->_states_count};
}

Complex const &
Regular
::echo() const
{
    return this->_states[0];
}

void
Regular
::apply_pulse(Quantity angle, Quantity phase)
{
    auto const T = operators::pulse(angle, phase);
    
    #pragma omp parallel for
    for(int order=0; order<this->_states_count; ++order)
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
Regular
::apply_time_interval(Quantity const & duration, Quantity const & gradient)
{
    // Note that since E does not depend on k, the E and S operators commute
    // and that E and D(k) also commute as they are diagonal matrices. The
    // only effect will be the relative order of D and S.
    // Since the diffusion operator relies on the "start" state k_1, we need
    // to apply the gradient operator after the diffusion operator. Otherwise
    // states would be dephased by D(k+Δk, Δk) instead of D(k, Δk)
    
    this->relaxation(duration);
    this->diffusion(duration, gradient);
    if(duration.magnitude != 0 && gradient.magnitude != 0)
    {
        this->shift();
    }
}

void
Regular
::shift()
{
    // TODO: resize factor
    if(3*this->_states_count >= this->_states.size())
    {
        this->_states.resize(this->_states.size()+3*100, 0);
    }
    
    // Shift positive F̃ states right
    for(int order=this->_states_count-1; order>=0; --order)
    {
        this->_states[3*(order+1)] = this->_states[3*order];
    }
    
    // Shift negative F̃^* states left
    for(int order=1; order<=this->_states_count; ++order)
    {
        this->_states[1+3*(order-1)] = this->_states[1+3*order];
    }
    
    // Update F̃_{+0} using F̃^*_{-0}
    this->_states[0] = std::conj(this->_states[1]);
    
    ++this->_states_count;
}

void
Regular
::relaxation(Quantity const & duration)
{
    if(this->species.get_R1().magnitude == 0 && this->species.get_R2().magnitude == 0)
    {
        return;
    }
    
    auto const E = operators::relaxation(this->species, duration);
    
    #pragma omp parallel for
    for(int order=0; order<this->_states_count; ++order)
    {
        this->_states[0+3*order] *= E.second;
        this->_states[1+3*order] *= E.second;
        this->_states[2+3*order] *= E.first;
    }
    
    this->_states[2] += 1.-E.first; // WARNING: assumes M0=1
}

void
Regular
::diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(this->species.get_D()[0].magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = sycomore::gamma*gradient*duration;
    
    #pragma omp parallel for
    for(int order=0; order<this->_states_count; ++order)
    {
        auto const k = order*delta_k;
        auto const D = operators::diffusion(this->species, duration, k, delta_k);
        this->_states[0+3*order] *= std::get<0>(D);
        this->_states[1+3*order] *= std::get<1>(D);
        this->_states[2+3*order] *= std::get<2>(D);
    }
}

}

}
