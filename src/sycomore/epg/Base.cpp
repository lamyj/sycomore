#include "Base.h"

#include <cmath>
#include <memory>

#include "sycomore/epg/operators.h"
#include "sycomore/epg/pool_model.h"
#include "sycomore/epg/pool_storage.h"
#include "sycomore/epg/simd_api.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"

namespace sycomore
{

namespace epg
{

Base
::Base(
    Species const & species, Magnetization const & initial_magnetization,
    unsigned int initial_size)
: _storage(new pool_storage::SinglePool(initial_size, 0)),
    _model(new pool_model::SinglePool(species, 0))
{
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_storage->F[0] = std::sqrt(2)*magnetization.p;
    this->_storage->F_star[0] = std::sqrt(2)*magnetization.m;
    this->_storage->Z[0] = magnetization.z;
    this->_model->M_z_eq = magnetization.z;
}

Base
::Base(Base const & other)
{
    this->operator=(other);
}

Base
::Base(Base && other)
: _storage(std::move(other._storage)), _model(std::move(other._model))
{
    // Nothing else.
}

Base &
Base
::operator=(Base const & other)
{
    if(std::dynamic_pointer_cast<pool_storage::SinglePool>(other._storage))
    {
        this->_storage = std::make_shared<pool_storage::SinglePool>(
            *std::dynamic_pointer_cast<pool_storage::SinglePool>(other._storage));
    }
    else if(std::dynamic_pointer_cast<pool_storage::Exchange>(other._storage))
    {
        this->_storage = std::make_shared<pool_storage::Exchange>(
            *std::dynamic_pointer_cast<pool_storage::Exchange>(other._storage));
    }
    else if(std::dynamic_pointer_cast<pool_storage::MagnetizationTransfer>(other._storage))
    {
        this->_storage = std::make_shared<pool_storage::MagnetizationTransfer>(
            *std::dynamic_pointer_cast<pool_storage::MagnetizationTransfer>(other._storage));
    }
    
    if(std::dynamic_pointer_cast<pool_model::SinglePool>(other._model))
    {
        this->_model = std::make_shared<pool_model::SinglePool>(
            *std::dynamic_pointer_cast<pool_model::SinglePool>(other._model));
    }
    else if(std::dynamic_pointer_cast<pool_model::Exchange>(other._model))
    {
        this->_model = std::make_shared<pool_model::Exchange>(
            *std::dynamic_pointer_cast<pool_model::Exchange>(other._model));
    }
    else if(std::dynamic_pointer_cast<pool_model::MagnetizationTransfer>(other._model))
    {
        this->_model = std::make_shared<pool_model::MagnetizationTransfer>(
            *std::dynamic_pointer_cast<pool_model::MagnetizationTransfer>(other._model));
    }
    
    return *this;
}

Base &
Base
::operator=(Base && other)
{
    this->_storage = std::move(other._storage);
    this->_model = std::move(other._model);
    
    return *this;
}

Species const &
Base
::get_species() const
{
    return this->_model->species;
}

void
Base
::set_species(Species const & species)
{
    this->_model->species = species;
}

std::vector<Complex>
Base
::state(std::size_t order) const
{
    return {
        this->_storage->F[order],
        this->_storage->F_star[order],
        this->_storage->Z[order]};
}

std::vector<Complex>
Base
::states() const
{
    std::vector<Complex> result(3*this->size());
    for(unsigned int order=0; order<this->size(); ++order)
    {
        result[3*order+0] = this->_storage->F[order];
        result[3*order+1] = this->_storage->F_star[order];
        result[3*order+2] = this->_storage->Z[order];
    }
    return result;
}

Complex const &
Base
::echo() const
{
    return this->_storage->F[0];
}

void
Base
::apply_pulse(Quantity angle, Quantity phase)
{
    simd_api::apply_pulse_single_pool(
        operators::pulse_single_pool(angle.magnitude, phase.magnitude), 
        *this->_storage, this->size());
}

void
Base
::relaxation(Quantity const & duration)
{
    if(
        this->_model->species.get_R1().magnitude == 0 
        && this->_model->species.get_R2().magnitude == 0)
    {
        return;
    }
    
    auto model = std::dynamic_pointer_cast<pool_model::SinglePool>(this->_model);
    auto const E = operators::relaxation_single_pool(
        model->species.get_R1().magnitude, model->species.get_R2().magnitude, 
        duration.magnitude);
    
    simd_api::relaxation_single_pool(
        E, *std::dynamic_pointer_cast<pool_storage::SinglePool>(this->_storage),
        this->size());
    
    this->_storage->Z[0] += this->_model->M_z_eq*(1.-E.first);
}

void
Base
::off_resonance(Quantity const & duration)
{
    auto const angle = 
        duration.magnitude * 2*M_PI*units::rad 
        * (
            this->delta_omega.magnitude
            + this->_model->species.get_delta_omega().magnitude);
    if(angle != 0)
    {
        auto const rotations = operators::phase_accumulation(angle);
        simd_api::off_resonance(
            rotations,
            this->_storage->F.data(),
            this->_storage->F_star.data(),
            this->_storage->Z.data(),
            this->size());
    }
}

}

}
