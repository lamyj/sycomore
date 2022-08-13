#ifndef _d635f223_7e0b_4f80_ab43_c03b499864f2
#define _d635f223_7e0b_4f80_ab43_c03b499864f2

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
protected:
    // NOTE: these need to be defined early so that we can bind the reference to
    // the "species" member.
    pool_storage::Base &_storage;
    pool_model::Base & _model;
public:
    Species & species;
    Real threshold=0;
    
    Base(
        Species const & species, Magnetization const & initial_magnetization,
        unsigned int initial_size);
    virtual ~Base();
};

}

}

#endif // _d635f223_7e0b_4f80_ab43_c03b499864f2
