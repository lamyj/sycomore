#include "pool_storage.h"

#include <utility>
#include <vector>
#include <xsimd/xsimd.hpp>
#include "sycomore/sycomore.h"

namespace sycomore
{

namespace epg
{

namespace pool_storage
{

Base
::Base(std::size_t size)
: F(size), F_star(size), Z(size)
{
    // Nothing else.
}

Base
::Base(std::size_t size, Complex const & value)
: F(size, value), F_star(size, value), Z(size, value)
{
    // Nothing else.
}


Exchange
::Exchange(std::size_t size)
: Base(size),
    F_a(this->F), F_star_a(this->F_star), Z_a(this->Z),
    F_b(size), F_star_b(size), Z_b(size)
{
    // Nothing else.
}

Exchange
::Exchange(std::size_t size, Complex const & value)
: Base(size, value),
    F_a(this->F), F_star_a(this->F_star), Z_a(this->Z),
    F_b(size, value), F_star_b(size, value), Z_b(size, value)
{
    // Nothing else.
}

Exchange
::Exchange(Exchange const & other)
: Base(other),
    F_a(this->F), F_star_a(this->F_star), Z_a(this->Z),
    F_b(other.F_b), F_star_b(other.F_star_b),
    Z_b(other.Z_b)
{
    // Nothing else.
}

Exchange
::Exchange(Exchange && other)
: Base(std::move(other)),
    F_a(this->F), F_star_a(this->F_star), Z_a(this->Z),
    F_b(std::move(other.F_b)), F_star_b(std::move(other.F_star_b)),
    Z_b(std::move(other.Z_b))
{
    // Nothing else.
}

Exchange &
Exchange
::operator=(Exchange const & other)
{
    this->Base::operator=(other);
    this->F_b = other.F_b;
    this->F_star_b = other.F_star_b;
    this->Z_b = other.Z_b;
    return *this;
}

Exchange &
Exchange
::operator=(Exchange && other)
{
    this->Base::operator=(std::move(other));
    this->F_b = std::move(other.F_b);
    this->F_star_b = std::move(other.F_star_b);
    this->Z_b = std::move(other.Z_b);
    return *this;
}

MagnetizationTransfer
::MagnetizationTransfer(std::size_t size)
: Base(size), Z_a(this->Z), Z_b(size)
{
    // Nothing else.
}

MagnetizationTransfer
::MagnetizationTransfer(std::size_t size, Complex const & value)
: Base(size), Z_a(this->Z), Z_b(size, value)
{
    // Nothing else.
}

MagnetizationTransfer
::MagnetizationTransfer(MagnetizationTransfer const & other)
: Base(other), Z_a(this->Z), Z_b(other.Z_b)
{
    // Nothing else.
}

MagnetizationTransfer
::MagnetizationTransfer(MagnetizationTransfer && other)
: Base(std::move(other)), Z_a(this->Z), Z_b(std::move(other.Z_b))
{
    // Nothing else.
}

MagnetizationTransfer &
MagnetizationTransfer
::operator=(MagnetizationTransfer const & other)
{
    this->Base::operator=(other);
    this->Z_b = other.Z_b;
    return *this;
}

MagnetizationTransfer &
MagnetizationTransfer
::operator=(MagnetizationTransfer && other)
{
    this->Base::operator=(std::move(other));
    this->Z_b = std::move(other.Z_b);
    return *this;
}

}

}

}
