#include "Base.h"

#include <cmath>

#include "sycomore/epg/pool_model.h"
#include "sycomore/epg/pool_storage.h"
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

}

}
