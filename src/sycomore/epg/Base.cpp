#include "Base.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "sycomore/epg/Model.h"
#include "sycomore/epg/operators.h"
#include "sycomore/epg/simd_api.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

Base
::Base(
    Species const & species, Magnetization const & M0,
    unsigned int initial_size)
: _model(species, M0, initial_size), _elapsed(0.)
{
    // Nothing else.
}

Base
::Base(
    Species const & species_a, Species const & species_b,
    Magnetization const & M0_a, Magnetization const & M0_b,
    Quantity const & k_a, Quantity const & delta_b,
    unsigned int initial_size)
: _model(species_a, species_b, M0_a, M0_b, k_a, delta_b, initial_size),
    _elapsed(0.)
{
    // Nothing else.
}

Base
::Base(
    Species const & species_a, Quantity const & R1_b_or_T1_b,
    Magnetization const & M0_a, Magnetization const & M0_b,
    Quantity const & k_a,
    unsigned int initial_size)
: _model(species_a, R1_b_or_T1_b, M0_a, M0_b, k_a, initial_size), _elapsed(0.)
{
    // Nothing else.
}

Model::Kind
Base
::kind() const
{
    return this->_model.kind;
}

std::size_t
Base
::pools() const
{
    return this->_model.pools;
}

Species const &
Base
::get_species(std::size_t pool) const
{
    return this->_model.species[pool];
}

void
Base
::set_species(std::size_t pool, Species const & species)
{
    this->_model.species[pool] = species;
}

void
Base
::set_species(Species const & species)
{
    this->set_species(0, species);
}

Real const &
Base
::get_M0(std::size_t pool) const
{
    return this->_model.M0[pool];
}

void
Base
::set_M0(std::size_t pool, Real const & M0)
{
    this->_model.M0[pool] = M0;
}

void
Base
::set_M0(Real const & M0)
{
    this->set_M0(0, M0);
}

Quantity const &
Base
::get_k(std::size_t pool) const
{
    return this->_model.k[pool];
}

void
Base
::set_k(std::size_t pool, Quantity const & k)
{
    this->_model.k[pool] = k;
}

Quantity const &
Base
::get_delta_b() const
{
    return this->_model.delta_b;
}

void
Base
::set_delta_b(Quantity const & delta_b)
{
    this->_model.delta_b = delta_b;
}

Base::State
Base
::state(std::size_t order) const
{
    Base::State result(3*this->_model.pools);
    for(std::size_t pool=0; pool < this->_model.pools; ++pool)
    {
        result[3*pool+0] = this->_model.F[pool][order];
        result[3*pool+1] = this->_model.F_star[pool][order];
        result[3*pool+2] = this->_model.Z[pool][order];
    }
    return result;
}

Base::States
Base
::states() const
{
    // We are keeping one set of orders for all pools: it is more logical to
    // order the dimensions as order > pool > (F, F*, Z)
    Base::States result(3*this->size()*this->_model.pools);
    auto it = result.begin();
    for(std::size_t order=0; order<this->size(); ++order)
    {
        for(std::size_t pool=0; pool < this->_model.pools; ++pool)
        {
            *(it++) = this->_model.F[pool][order];
            *(it++) = this->_model.F_star[pool][order];
            *(it++) = this->_model.Z[pool][order];
        }
    }
    
    return result;
}

Quantity
Base
::elapsed() const
{
    return this->_elapsed*sycomore::units::s;
}

Complex const &
Base
::echo(std::size_t pool) const
{
    return this->_model.F[pool][0];
}

void
Base
::apply_pulse(Quantity const & angle, Quantity const & phase)
{
    if(this->_model.kind == Model::SinglePool)
    {
        auto const T = operators::pulse_single_pool(
        angle.magnitude, phase.magnitude);
        simd_api::apply_pulse_single_pool(T, this->_model, this->size());
    }
    else if(this->_model.kind == Model::Exchange)
    {
        auto const T = operators::pulse_exchange(
            angle.magnitude, phase.magnitude,
            angle.magnitude, phase.magnitude);
        simd_api::apply_pulse_exchange(T, this->_model, this->size());
    }
    else
    {
        throw std::runtime_error("Invalid model");
    }
}

void
Base
::apply_pulse(
    Quantity const & angle_a, Quantity const & phase_a,
    Quantity const & angle_b, Quantity const & phase_b)
{
    if(this->_model.kind != Model::Exchange)
    {
        throw std::runtime_error("Invalid model");
    }
    
    auto const T = operators::pulse_exchange(
        angle_a.magnitude, phase_a.magnitude,
        angle_b.magnitude, phase_b.magnitude);
    simd_api::apply_pulse_exchange(T, this->_model, this->size());
}

void
Base
::apply_pulse(
    Quantity const & angle_a, Quantity const & phase_a, Real saturation)
{
    if(this->_model.kind != Model::MagnetizationTransfer)
    {
        throw std::runtime_error("Invalid model");
    }
    
    auto const T = operators::pulse_magnetization_transfer(
        angle_a.magnitude, phase_a.magnitude, saturation);
    simd_api::apply_pulse_magnetization_transfer(T, this->_model, this->size());
}

void
Base
::relaxation(Quantity const & duration)
{
    if(this->_model.kind == Model::SinglePool)
    {
        auto const & species = this->_model.species[0];
        if(species.get_R1().magnitude == 0 && species.get_R2().magnitude == 0)
        {
            return;
        }
        auto const E = operators::relaxation_single_pool(
            species.get_R1().magnitude, species.get_R2().magnitude, 
            duration.magnitude);
        simd_api::relaxation_single_pool(E, this->_model, this->size());
        this->_model.Z[0][0] += this->_model.M0[0]*(1.-E.first);
    }
    else if(this->_model.kind == Model::Exchange)
    {
        auto const E = operators::relaxation_exchange(
            this->_model.species[0].get_R1().magnitude,
            this->_model.species[0].get_R2().magnitude,
            this->_model.species[1].get_R1().magnitude,
            this->_model.species[1].get_R2().magnitude,
            this->_model.k[0].magnitude, this->_model.k[1].magnitude,
            this->_model.delta_b.magnitude,
            this->_model.M0[0], this->_model.M0[1],
            duration.magnitude);
        simd_api::relaxation_exchange(
            std::get<0>(E), std::get<1>(E), this->_model, this->size());
        
        auto const & recovery = std::get<2>(E);
        this->_model.Z[0][0] += recovery[0];
        this->_model.Z[1][0] += recovery[1];
    }
    else if(this->_model.kind == Model::MagnetizationTransfer)
    {
        auto const E = operators::relaxation_magnetization_transfer(
            this->_model.species[0].get_R1().magnitude,
            this->_model.species[0].get_R2().magnitude,
            this->_model.species[1].get_R1().magnitude,
            this->_model.k[0].magnitude, this->_model.k[1].magnitude,
            this->_model.M0[0], this->_model.M0[1],
            duration.magnitude);
        simd_api::relaxation_magnetization_transfer(
            std::get<0>(E), std::get<1>(E), this->_model, this->size());
        
        auto const & recovery = std::get<2>(E);
        this->_model.Z[0][0] += recovery[0];
        this->_model.Z[1][0] += recovery[1];
    }
    else
    {
        throw std::runtime_error("Invalid model");
    }
}

void
Base
::off_resonance(Quantity const & duration)
{
    for(std::size_t pool=0; pool<this->_model.pools; ++pool)
    {
        auto const angle = 
            duration.magnitude * 2*M_PI*units::rad 
            * (
                this->delta_omega.magnitude
                + this->_model.species[pool].get_delta_omega().magnitude);
        if(angle != 0)
        {
            auto const Omega = operators::phase_accumulation(angle);
            simd_api::off_resonance(
                Omega, this->_model.F[pool], this->_model.F_star[pool],
                this->size());
        }
    }
}

}

}
