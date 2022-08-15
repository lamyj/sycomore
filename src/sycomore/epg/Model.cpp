#include "Model.h"

#include <stdexcept>
#include <vector>
#include <utility>

#include <xsimd/xsimd.hpp>

#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"

namespace sycomore
{

namespace epg
{

Model
::Model(Species const & species, Magnetization M0, std::size_t initial_size)
: kind(SinglePool), pools(1),
    species({species}), M0(pools), k(0), delta_b(0*units::Hz),
    F(pools), F_star(pools), Z(pools)
{
    this->_initialize({M0}, initial_size);
}

Model
::Model(
    Species const & species_a, Species const & species_b,
    Magnetization const & M0_a, Magnetization const & M0_b,
    Quantity const & k_a, Quantity const & delta_b,
    std::size_t initial_size)
: kind(Exchange), pools(2),
    species({species_a, species_b}), M0(pools), k(pools), delta_b(delta_b),
    F(pools), F_star(pools), Z(pools)
{
    this->_initialize({M0_a, M0_b}, initial_size);
    this->k = {k_a, k_a*this->M0[0]/this->M0[1]};
}

Model
::Model(
    Species const & species_a, Quantity const & R1_b_or_T1_b,
    Magnetization const & M0_a, Magnetization const & M0_b,
    Quantity const & k_a,
    std::size_t initial_size)
: kind(MagnetizationTransfer), pools(2),
    species({species_a, Species(R1_b_or_T1_b, 0*units::Hz)}),
    M0(pools), k(pools), delta_b(0*units::Hz),
    F(pools), F_star(pools), Z(pools)
{
    this->_initialize({M0_a, M0_b}, initial_size);
    this->k = {k_a, k_a*this->M0[0]/this->M0[1]};
}

Model &
Model
::operator=(Model const & other)
{
    if(this->kind != other.kind)
    {
        throw std::runtime_error("Invalid model assignment");
    }
    this->species = other.species;
    this->M0 = other.M0;
    this->k = other.k;
    this->delta_b = other.delta_b;
    
    this->F = other.F;
    this->F_star = other.F_star;
    this->Z = other.Z;
    
    return *this;
}

Model &
Model
::operator=(Model && other)
{
    if(this->kind != other.kind)
    {
        throw std::runtime_error("Invalid model assignment");
    }
    this->species = std::move(other.species);
    this->M0 = std::move(other.M0);
    this->k = std::move(other.k);
    this->delta_b = std::move(other.delta_b);
    
    this->F = std::move(other.F);
    this->F_star = std::move(other.F_star);
    this->Z = std::move(other.Z);
    
    return *this;
}

void
Model
::_initialize(std::vector<Magnetization> const & M0, std::size_t initial_size)
{
    for(auto & item: this->F)
    {
        item.resize(initial_size);
    }
    for(auto & item: this->F_star)
    {
        item.resize(initial_size);
    }
    for(auto & item: this->Z)
    {
        item.resize(initial_size);
    }
    
    for(std::size_t i=0; i<this->pools; ++i)
    {
        auto const M0_c = as_complex_magnetization(M0[i]);
        this->F[i][0] = std::sqrt(2)*M0_c.p;
        this->F_star[i][0] = std::sqrt(2)*M0_c.m;
        this->Z[i][0] = M0_c.z;
        this->M0[i] = M0_c.z;
    }
}

}

}
