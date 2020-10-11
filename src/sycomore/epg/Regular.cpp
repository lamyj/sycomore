#include "Regular.h"

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
    unsigned int initial_size, 
    Quantity const & unit_gradient_area, double gradient_tolerance)
: species(species),
    _F(initial_size, 0), _F_star(initial_size, 0), _Z(initial_size, 0),
    _unit_gradient_area(unit_gradient_area), _gradient_tolerance(gradient_tolerance)
{
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_F[0] = std::sqrt(2)*magnetization.p;
    this->_F_star[0] = std::sqrt(2)*magnetization.m;
    this->_Z[0] = magnetization.z;
    
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
    return {this->_F[order], this->_F_star[order], this->_Z[order]};
}

std::vector<Complex>
Regular
::states() const
{
    std::vector<Complex> result(3*this->_states_count);
    for(unsigned int order=0; order<this->_states_count; ++order)
    {
        result[3*order+0] = this->_F[order];
        result[3*order+1] = this->_F_star[order];
        result[3*order+2] = this->_Z[order];
    }
    return result;
}

Complex const &
Regular
::echo() const
{
    return this->_F[0];
}

void
Regular
::apply_pulse(Quantity angle, Quantity phase)
{
    auto const T = operators::pulse(angle.magnitude, phase.magnitude);
    
    Complex F, F_star, Z;
    
    for(int order=0; order<this->_states_count; ++order)
    {
        F = 
            T[3*0+0] * this->_F[order]
            + T[3*0+1] * this->_F_star[order]
            + T[3*0+2] * this->_Z[order];
        F_star = 
            T[3*1+0] * this->_F[order]
            + T[3*1+1] * this->_F_star[order]
            + T[3*1+2] * this->_Z[order];
        Z = 
            T[3*2+0] * this->_F[order]
            + T[3*2+1] * this->_F_star[order]
            + T[3*2+2] * this->_Z[order];
        
        this->_F[order] = F;
        this->_F_star[order] = F_star;
        this->_Z[order] = Z;
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
    this->bulk_motion(duration, gradient);
    if(duration.magnitude != 0 && gradient.magnitude != 0)
    {
        if(this->_unit_gradient_area.magnitude != 0)
        {
            this->shift(duration, gradient);
        }
        else
        {
            this->shift();
        }
    }
    if(
        duration.magnitude != 0 
        && (this->delta_omega.magnitude != 0 || species.get_delta_omega().magnitude != 0))
    {
        this->off_resonance(duration);
    }
}

void
Regular
::apply_time_interval(TimeInterval const & interval)
{
    this->apply_time_interval(
        interval.get_duration(), interval.get_gradient_amplitude()[0]);
}

void
Regular
::shift()
{
    this->_shift(1);
}

void
Regular
::shift(Quantity const & duration, Quantity const & gradient)
{
    auto const area = duration*gradient;
    auto const epsilon = this->_gradient_tolerance*this->_unit_gradient_area.magnitude;
    auto const remainder = std::remainder(
        area.magnitude, this->_unit_gradient_area.magnitude);
    
    if(std::abs(remainder) >= epsilon)
    {
        throw std::runtime_error(
            "Gradient is not a interger multiple of unit gradient");
    }
    
    int n = std::round(area/this->_unit_gradient_area);
    
    this->_shift(n);
}

void
Regular
::relaxation(Quantity const & duration)
{
    if(this->species.get_R1().magnitude == 0 && this->species.get_R2().magnitude == 0)
    {
        return;
    }
    
    auto const E = operators::relaxation(
        this->species.get_R1().magnitude, this->species.get_R2().magnitude, 
        duration.magnitude);
    
    for(int order=0; order<this->_states_count; ++order)
    {
        this->_F[order] *= E.second;
        this->_F_star[order] *= E.second;
        this->_Z[order] *= E.first;
    }
    
    this->_Z[0] += 1.-E.first; // WARNING: assumes M0=1
}

void
Regular
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
    
    #pragma omp parallel for schedule(static)
    for(int order=0; order<this->_states_count; ++order)
    {
        auto const k = order*delta_k;
        auto const D = operators::diffusion(
            this->species.get_D()[0].magnitude, duration.magnitude, k, delta_k);
        this->_F[order] *= std::get<0>(D);
        this->_F_star[order] *= std::get<1>(D);
        this->_Z[order] *= std::get<2>(D);
    }
}

void
Regular
::off_resonance(Quantity const & duration)
{
    auto const angle = 
        duration * 2*M_PI*units::rad 
        * (this->delta_omega+this->species.get_delta_omega());
    if(angle.magnitude != 0)
    {
        auto const rotations = operators::phase_accumulation(angle.magnitude);
        
        for(int order=0; order<this->_states_count; ++order)
        {
            this->_F[order] *= rotations.first;
            this->_F_star[order] *= rotations.second;
            // Z̃ states are unaffected
        }
    }
}

void
Regular
::bulk_motion(Quantity const & duration, Quantity const & gradient)
{
    if(this->velocity.magnitude == 0)
    {
        return;
    }
    
    auto const delta_k = (sycomore::gamma*gradient*duration).magnitude;
    if(delta_k == 0)
    {
        return;
    }
    
    #pragma omp parallel for schedule(static)
    for(int order=0; order<this->_states_count; ++order)
    {
        auto const k = order*delta_k;
        auto const J = operators::bulk_motion(
            this->velocity.magnitude, duration.magnitude, k, delta_k);
        this->_F[order] *= std::get<0>(J);
        this->_F_star[order] *= std::get<1>(J);
        this->_Z[order] *= std::get<2>(J);
    }
}

Quantity const &
Regular
::unit_gradient_area() const
{
    return this->_unit_gradient_area;
}

double
Regular
::gradient_tolerance() const
{
    return this->_gradient_tolerance;
}

void
Regular
::_shift(int n)
{
    if(std::abs(n)>1)
    {
        for(int i=0; i<std::abs(n); ++i)
        {
            this->_shift(n>0?+1:-1);
        }
    }
    else if(n==1 || n==-1)
    {
        if(this->_states_count >= this->_F.size())
        {
            this->_F.resize(this->_F.size()+100, 0);
            this->_F_star.resize(this->_F_star.size()+100, 0);
            this->_Z.resize(this->_Z.size()+100, 0);
        }
        
        if(n == +1)
        {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    // Shift positive F̃ states right
                    std::copy_backward(
                        this->_F.begin(), this->_F.begin()+this->_states_count, 
                        this->_F.begin()+this->_states_count+1);
                }
                #pragma omp section
                {
                    // Shift negative F̃^* states left
                    std::copy(
                        this->_F_star.begin()+1, 
                        this->_F_star.begin()+this->_states_count+1, 
                        this->_F_star.begin());
                }
            }
            
            // Update F̃_{+0} using F̃^*_{-0}
            this->_F[0] = std::conj(this->_F_star[0]);
        }
        else
        {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    // Shift negative F̃^* states right
                    std::copy_backward(
                        this->_F_star.begin(), 
                        this->_F_star.begin()+this->_states_count, 
                        this->_F_star.begin()+this->_states_count+1);
                }
                #pragma omp section
                {
                    // Shift positive F̃ states left
                    std::copy(
                        this->_F.begin()+1, 
                        this->_F.begin()+this->_states_count+1, 
                        this->_F.begin());
                }
            }
            
            // Update F̃^*_{-0} using F̃_{+0}
            this->_F_star[0] = std::conj(this->_F[0]);
        }
        
        ++this->_states_count;
        
        // Remove empty states with high order.
        while(
            this->_F[this->_states_count-1] == 0. 
            && this->_F_star[this->_states_count-1] == 0.
            && this->_Z[this->_states_count-1] == 0.)
        {
            --this->_states_count;
        }
    }
    else
    { 
        // n==0, nothing to do.
    }
}

}

}
