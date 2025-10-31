#ifndef _f4271617_4909_4e41_bd07_ee1eb44f75ab
#define _f4271617_4909_4e41_bd07_ee1eb44f75ab

#include <stdexcept>

#include "QuantityInterface.h"

namespace sycomore
{

template<typename TDerived>
void
QuantityInterface<TDerived>
::check_dimensions(Dimensions const & dimensions, std::string const & message) const
{
    if(this->derived_cast().dimensions != dimensions)
    {
        throw std::runtime_error(message);
    }
}

template<typename TDerived>
template<typename T>
void
QuantityInterface<TDerived>
::check_dimensions(T const & other, std::string const & message) const
{
    this->check_dimensions(other.dimensions);
}

template<typename TDerived>
TDerived &
QuantityInterface<TDerived>
::derived_cast()
{
    return static_cast<TDerived &>(*this);
}

template<typename TDerived>
TDerived const &
QuantityInterface<TDerived>
::derived_cast() const
{
    return static_cast<TDerived const &>(*this);
}

template<typename TDerived>
bool
QuantityInterface<TDerived>
::operator==(TDerived const & right) const
{
    auto const & left = this->derived_cast();
    return
        left.magnitude == right.magnitude
        && left.dimensions == right.dimensions;
}

template<typename TDerived>
bool
QuantityInterface<TDerived>
::operator!=(TDerived const & right) const
{
    auto const & left = this->derived_cast();
    return
        left.magnitude != right.magnitude
        || left.dimensions != right.dimensions;
}

template<typename TDerived>
template<typename T>
TDerived &
QuantityInterface<TDerived>
::operator+=(T const & right)
{
    this->check_dimensions(right, "Addition requires same dimensions");
    auto & left = this->derived_cast();
    left.magnitude += right.magnitude;
    return left;
}

template<typename TDerived>
template<typename T>
TDerived &
QuantityInterface<TDerived>
::operator-=(T const & right)
{
    this->check_dimensions(right, "Subtraction requires same dimensions");
    auto & left = this->derived_cast();
    left.magnitude -= right.magnitude;
    return left;
}

template<typename TDerived>
template<typename T>
TDerived &
QuantityInterface<TDerived>
::fill(T const & x)
{
    auto & self = this->derived_cast();
    self.magnitude.fill(x.magnitude);
    self.dimensions = x.dimensions;
    return self;
}

template<typename TDerived>
std::ostream & operator<<(
    std::ostream & stream, QuantityInterface<TDerived> const & q)
{
    auto const & q_ = q.derived_cast();
    stream << q_.magnitude << " " << q_.dimensions << std::endl;
    return stream;
}

}

#endif // _f4271617_4909_4e41_bd07_ee1eb44f75ab
