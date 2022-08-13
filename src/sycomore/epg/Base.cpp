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
: _storage(*new pool_storage::SinglePool(initial_size, 0)),
    _model(*new pool_model::SinglePool(species, 0)), species(_model.species)
{
    auto const magnetization = as_complex_magnetization(initial_magnetization);
    this->_storage.F[0] = std::sqrt(2)*magnetization.p;
    this->_storage.F_star[0] = std::sqrt(2)*magnetization.m;
    this->_storage.Z[0] = magnetization.z;
    this->_model.M_z_eq = magnetization.z;
}

Base
::~Base()
{
    delete &this->_storage;
    delete &this->_model;
}

}

}
