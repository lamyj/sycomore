#ifndef _961e2ca9_d932_4fa0_b947_437f66689738
#define _961e2ca9_d932_4fa0_b947_437f66689738

#include <sstream>
#include <stdexcept>

#include "sycomore/Quantity.h"
#include "sycomore/Species.h"

namespace sycomore
{

namespace epg
{

/**
 * @brief Model parameters for single-pool, two pools with exchange, and two
 * pools with magnetization transfer.
 */
namespace pool_model
{

/// @brief Single pool: one species and one equilibrium Z magnetization.
struct Base
{
    Species species;
    Real M_z_eq;
    
    Base(Species const & species, Real M_z_eq);
    
    Base(Base const &) = default;
    Base(Base &&) = default;
    Base & operator=(Base const &) = default;
    Base & operator=(Base &&) = default;
    virtual ~Base() = default;
};

using SinglePool = Base;

/**
 * @brief Two pools with exchange: two separate species, equilibrium 
 * magnetizations and exchange rates, and frequency offset.
 */
struct Exchange: public Base
{
    Species & species_a;
    Real & M_z_eq_a;
    
    Species species_b;
    Real M_z_eq_b;
    Quantity k_a, k_b;
    Quantity delta_b;

    Exchange(
        Species const & species_a, Species const & species_b,
        Real M_z_eq_a, Real M_z_eq_b,
        Quantity const & k_a, Quantity const & delta_b);
    
    Exchange(Exchange const & other);
    Exchange(Exchange && other);
    Exchange & operator=(Exchange const & other);
    Exchange & operator=(Exchange && other);
    
    virtual ~Exchange() = default;
};

/**
 * @brief Two pools with magnetization transfer: the second pool is assumed to
 * have no transversal magnetization.
 */
struct MagnetizationTransfer: public Base
{
    Species & species_a;
    Real & M_z_eq_a;
    
    Quantity R1_b;
    Real M_z_eq_b;
    Quantity k_a, k_b;

    MagnetizationTransfer(
        Species const & a, Quantity const & R1_b_or_T1_b,
        Real M_z_eq_a, Real M_z_eq_b,
        Quantity const & k_a);
    
    MagnetizationTransfer(MagnetizationTransfer const & other);
    MagnetizationTransfer(MagnetizationTransfer && other);
    MagnetizationTransfer & operator=(MagnetizationTransfer const & other);
    MagnetizationTransfer & operator=(MagnetizationTransfer && other);
    
    virtual ~MagnetizationTransfer() = default;
};

}

}

}

#endif // _961e2ca9_d932_4fa0_b947_437f66689738
