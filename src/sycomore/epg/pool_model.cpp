#include "pool_model.h"

#include <sstream>
#include <stdexcept>

#include "sycomore/Quantity.h"
#include "sycomore/Species.h"

namespace sycomore
{

namespace epg
{

namespace pool_model
{

Base
::Base(Species const & species, Real M_z_eq)
: species(species), M_z_eq(M_z_eq)
{
    // Nothing else.
}

Exchange
::Exchange(
    Species const & species_a, Species const & species_b,
    Real M_z_eq_a, Real M_z_eq_b,
    Quantity const & k_a, Quantity const & delta_b)
: Base(species_a, M_z_eq_a), 
    species_a(this->species), species_b(species_b),
    M_z_eq_a(this->M_z_eq), M_z_eq_b(M_z_eq_b),
    k_a(k_a), k_b(k_a*M_z_eq_a/M_z_eq_b), delta_b(delta_b)
{
    // Nothing else.
}

Exchange
::Exchange(Exchange const & other)
: Base(other), 
    species_a(this->species), species_b(other.species_b),
    M_z_eq_a(this->M_z_eq), M_z_eq_b(other.M_z_eq_b),
    k_a(other.k_a), k_b(this->k_a*this->M_z_eq_a/this->M_z_eq_b),
    delta_b(other.delta_b)
{
    // Nothing else.
}

Exchange
::Exchange(Exchange && other)
: Base(std::move(other)),
    species_a(this->species), species_b(std::move(other.species_b)),
    M_z_eq_a(this->M_z_eq), M_z_eq_b(std::move(other.M_z_eq_b)),
    k_a(std::move(other.k_a)), k_b(this->k_a*this->M_z_eq_a/this->M_z_eq_b),
    delta_b(std::move(other.delta_b))
{
    // Nothing else.
}

Exchange &
Exchange
::operator=(Exchange const & other)
{
    this->Base::operator=(other);
    this->species_b = other.species_b;
    this->M_z_eq_b = other.M_z_eq_b;
    this->k_a = other.k_a;
    this->k_b = this->k_a*this->M_z_eq_a/this->M_z_eq_b;
    this->delta_b = other.delta_b;
    
    return *this;
}

Exchange &
Exchange
::operator=(Exchange && other)
{
    this->Base::operator=(std::move(other));
    this->species_b = std::move(other.species_b);
    this->M_z_eq_b = std::move(other.M_z_eq_b);
    this->k_a = std::move(other.k_a);
    this->k_b = this->k_a*this->M_z_eq_a/this->M_z_eq_b;
    this->delta_b = std::move(other.delta_b);
    
    return *this;
}

MagnetizationTransfer
::MagnetizationTransfer(
    Species const & species_a, Quantity const & R1_b_or_T1_b,
    Real M_z_eq_a, Real M_z_eq_b,
    Quantity const & k_a)
: Base(species_a, M_z_eq_a), 
    species_a(this->species), 
    M_z_eq_a(this->M_z_eq), M_z_eq_b(M_z_eq_b),
    k_a(k_a), k_b(k_a*M_z_eq_a/M_z_eq_b)
{
    if(R1_b_or_T1_b.dimensions == Frequency)
    {
        this->R1_b = R1_b_or_T1_b;
    }
    else if(R1_b_or_T1_b.dimensions == Time)
    {
        this->R1_b = 1/R1_b_or_T1_b;
    }
    else
    {
        std::ostringstream message;
        message
            << "R1b must be duration or frequency, not "
            << R1_b_or_T1_b.dimensions;
        throw std::runtime_error(message.str());
    }
}

MagnetizationTransfer
::MagnetizationTransfer(MagnetizationTransfer const & other)
: Base(other),
    species_a(this->species), M_z_eq_a(this->M_z_eq),
    R1_b(other.R1_b), M_z_eq_b(other.M_z_eq_b),
    k_a(other.k_a), k_b(other.k_b)
{
    // Nothing else.
}

MagnetizationTransfer
::MagnetizationTransfer(MagnetizationTransfer && other)
: Base(std::move(other)),
    species_a(this->species), M_z_eq_a(this->M_z_eq),
    R1_b(std::move(other.R1_b)), M_z_eq_b(std::move(other.M_z_eq_b)),
    k_a(std::move(other.k_a)), k_b(std::move(other.k_b))
{
    // Nothing else.
}

MagnetizationTransfer &
MagnetizationTransfer
::operator=(MagnetizationTransfer const & other)
{
    this->Base::operator=(other);
    this->R1_b = other.R1_b;
    this->M_z_eq_b = other.M_z_eq_b;
    this->k_a = other.k_a;
    this->k_b = other.k_b;
    
    return *this;
}

MagnetizationTransfer &
MagnetizationTransfer
::operator=(MagnetizationTransfer && other)
{
    this->Base::operator=(std::move(other));
    this->R1_b = std::move(other.R1_b);
    this->M_z_eq_b = std::move(other.M_z_eq_b);
    this->k_a = std::move(other.k_a);
    this->k_b = std::move(other.k_b);
    
    return *this;
}


}

}

}
