#include "Regular.h"

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
    unsigned int initial_size, Quantity gamma)
: species(species), gamma(gamma), _magnetization(3*initial_size, 0)
{
    // Store magnetization as lines of F̃_k, F̃^*_{-k}, Z̃_k
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_magnetization[0] = std::sqrt(2)*magnetization.p;
    this->_magnetization[1] = std::sqrt(2)*magnetization.m;
    this->_magnetization[2] = magnetization.z;
    
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
::magnetization(std::size_t state) const
{
    return {
        this->_magnetization[3*state], 
        this->_magnetization[3*state+1], 
        this->_magnetization[3*state+2]};
}

Complex const &
Regular
::echo() const
{
    return this->_magnetization[0];
}

void
Regular
::apply_pulse(Quantity angle, Quantity phase)
{
    auto const T = operators::pulse(angle, phase);
    
    #pragma omp parallel for
    for(int s=0; s<this->_states_count; ++s)
    {
        std::vector<Complex> result(3, 0);
        for(int r=0; r<3; ++r)
        {
            for(int c=0; c<3; ++c)
            {
                result[r] += T[3*r+c] * this->_magnetization[3*s+c];
            }
        }
        
        std::memcpy(
            this->_magnetization.data()+3*s, result.data(), 3*sizeof(Complex));
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
    
    this->apply_relaxation(duration);
    this->apply_diffusion(duration, gradient);
    if(duration.magnitude != 0 && gradient.magnitude != 0)
    {
        this->apply_gradient();
    }
}

void
Regular
::apply_gradient()
{
    // TODO: resize factor
    if(3*this->_states_count >= this->_magnetization.size())
    {
        this->_magnetization.resize(this->_magnetization.size()+3*100, 0);
    }
    
    // Shift positive F̃ states right
    for(int s=this->_states_count-1; s>=0; --s)
    {
        this->_magnetization[3*(s+1)] = this->_magnetization[3*s];
    }
    
    // Shift negative F̃^* states left
    for(int s=1; s<=this->_states_count; ++s)
    {
        this->_magnetization[1+3*(s-1)] = this->_magnetization[1+3*s];
    }
    
    // Update F̃_{+0} using F̃^*_{-0}
    this->_magnetization[0] = std::conj(this->_magnetization[1]);
    
    ++this->_states_count;
}

void
Regular
::apply_relaxation(Quantity const & duration)
{
    if(this->species.get_R1().magnitude == 0 && this->species.get_R2().magnitude == 0)
    {
        return;
    }
    
    auto const E = operators::relaxation(this->species, duration);
    
    #pragma omp parallel for
    for(int s=0; s<this->_states_count; ++s)
    {
        this->_magnetization[0+3*s] *= E.second;
        this->_magnetization[1+3*s] *= E.second;
        this->_magnetization[2+3*s] *= E.first;
    }
    
    this->_magnetization[2] += 1.-E.first; // WARNING: assumes M0=1
}

void
Regular
::apply_diffusion(Quantity const & duration, Quantity const & gradient)
{
    if(this->species.get_D().magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = this->gamma*gradient*duration;
    
    #pragma omp parallel for
    for(int s=0; s<this->_states_count; ++s)
    {
        auto const k = s*delta_k;
        auto const D = operators::diffusion(this->species, duration, k, delta_k);
        this->_magnetization[0+3*s] *= std::get<0>(D);
        this->_magnetization[1+3*s] *= std::get<1>(D);
        this->_magnetization[2+3*s] *= std::get<2>(D);
    }
}

}

}
