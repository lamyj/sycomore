#ifndef _d635f223_7e0b_4f80_ab43_c03b499864f2
#define _d635f223_7e0b_4f80_ab43_c03b499864f2

#include <memory>

#include "sycomore/epg/pool_model.h"
#include "sycomore/epg/pool_storage.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"

namespace sycomore
{

namespace epg
{

/// @brief Base class for all EPG models.
class Base
{
public:
    Real threshold=0;
    
    Base(
        Species const & species, Magnetization const & initial_magnetization,
        unsigned int initial_size);
    Base(Base const & other);
    Base(Base && other);
    Base & operator=(Base const & other);
    Base & operator=(Base && other);
    virtual ~Base() = default;
    
    Species const & get_species() const;
    void set_species(Species const & species);
protected:
    std::shared_ptr<pool_storage::Base> _storage;
    std::shared_ptr<pool_model::Base> _model;
};

}

}

#endif // _d635f223_7e0b_4f80_ab43_c03b499864f2
