#ifndef _f4271617_4909_4e41_bd07_ee1eb44f75ab
#define _f4271617_4909_4e41_bd07_ee1eb44f75ab

#include <cmath>
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
    return !this->operator==(right);
}

template<typename TDerived>
template<typename T, enable_if_quantity<T>>
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
template<typename T, enable_if_not_quantity<T> >
TDerived &
QuantityInterface<TDerived>
::operator+=(T const & right)
{
    this->check_dimensions(Dimensionless, "Addition requires same dimensions");
    auto & left = this->derived_cast();
    left.magnitude += right;
    return left;
}

template<typename TDerived>
template<typename T, enable_if_quantity<T>>
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
template<typename T, enable_if_not_quantity<T> >
TDerived &
QuantityInterface<TDerived>
::operator-=(T const & right)
{
    this->check_dimensions(Dimensionless, "Subtraction requires same dimensions");
    auto & left = this->derived_cast();
    left.magnitude -= right;
    return left;
}

template<typename TDerived>
template<typename T, enable_if_quantity<T>>
TDerived &
QuantityInterface<TDerived>
::operator*=(T const & right)
{
    auto & left = this->derived_cast();
    left.magnitude *= right.magnitude;
    left.dimensions *= right.dimensions;
    return left;
}

template<typename TDerived>
template<typename T, enable_if_not_quantity<T>>
TDerived &
QuantityInterface<TDerived>
::operator*=(T const & right)
{
    auto & left = this->derived_cast();
    left.magnitude *= right;
    return left;
}

template<typename TDerived>
template<typename T, enable_if_quantity<T>>
TDerived &
QuantityInterface<TDerived>
::operator/=(T const & right)
{
    auto & left = this->derived_cast();
    left.magnitude /= right.magnitude;
    left.dimensions /= right.dimensions;
    return left;
}

template<typename TDerived>
template<typename T, enable_if_not_quantity<T>>
TDerived &
QuantityInterface<TDerived>
::operator/=(T const & right)
{
    auto & left = this->derived_cast();
    left.magnitude /= right;
    return left;
}

template<typename TDerived>
template<typename T, enable_if_quantity<T>>
TDerived &
QuantityInterface<TDerived>
::operator%=(T const & right)
{
    auto & left = this->derived_cast();
    if(right.dimensions == Dimensionless || left.dimensions == right.dimensions)
    {
        return this->operator%=(right.magnitude);
    }
    else
    {
        throw std::runtime_error("Modulo requires same dimensions");
    }
    return left;
}

template<typename TDerived>
template<typename T, enable_if_not_quantity<T>>
TDerived &
QuantityInterface<TDerived>
::operator%=(T const & right)
{
    auto & left = this->derived_cast();
    using std::fmod;
    using xt::fmod;
    left.magnitude = fmod(left.magnitude, right);
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

template<typename T, enable_if_quantity<T>>
bool operator==(typename T::Container const & left, T const & right)
{
    return right.operator==(left);
}

template<typename T, enable_if_quantity<T>>
bool operator!=(typename T::Container const & left, T const & right)
{
    return !(left == right);
}

template<typename T, enable_if_quantity<T>>
OwningType<T> operator+(T x)
{
    return {+x.magnitude, x.dimensions};
}

template<typename T, enable_if_quantity<T>>
OwningType<T> operator-(T const & x)
{
    return {-x.magnitude, x.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<T1, T2> operator+(T1 const & left, T2 const & right)
{
    left.check_dimensions(right, "Addition requires same dimensions");
    return {left.magnitude + right.magnitude, left.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_not_quantity<T2>>
CommonQuantityType<T1, QuantityContainer<T2>>
operator+(T1 const & left, T2 const & right)
{
    left.check_dimensions(Dimensionless, "Addition requires same dimensions");
    return {left.magnitude + right, left.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_not_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<QuantityContainer<T1>, T2>
operator+(T1 const & left, T2 const & right)
{
    right.check_dimensions(Dimensionless, "Addition requires same dimensions");
    return {left + right.magnitude, right.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<T1, T2> operator-(T1 const & left, T2 const & right)
{
    left.check_dimensions(right, "Subtraction requires same dimensions");
    return {left.magnitude - right.magnitude, left.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_not_quantity<T2>>
CommonQuantityType<T1, QuantityContainer<T2>>
operator-(T1 const & left, T2 const & right)
{
    left.check_dimensions(Dimensionless, "Subtraction requires same dimensions");
    return {left.magnitude - right, left.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_not_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<QuantityContainer<T1>, T2>
operator-(T1 const & left, T2 const & right)
{
    right.check_dimensions(Dimensionless, "Subtraction requires same dimensions");
    return {left - right.magnitude, right.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<T1, T2> operator*(T1 const & left, T2 const & right)
{
    return {left.magnitude * right.magnitude, left.dimensions * right.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_not_quantity<T2>>
CommonQuantityType<T1, QuantityContainer<T2>>
operator*(T1 const & left, T2 const & right)
{
    return {left.magnitude * right, left.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_not_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<QuantityContainer<T1>, T2>
operator*(T1 const & left, T2 const & right)
{
    return {left * right.magnitude, right.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<T1, T2> operator/(T1 const & left, T2 const & right)
{
    return {left.magnitude / right.magnitude, left.dimensions / right.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_not_quantity<T2>>
CommonQuantityType<T1, QuantityContainer<T2>>
operator/(T1 const & left, T2 const & right)
{
    return {left.magnitude / right, left.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_not_quantity<T1>, enable_if_quantity<T2>>
CommonQuantityType<QuantityContainer<T1>, T2>
operator/(T1 const & left, T2 const & right)
{
    return {left / right.magnitude, std::pow(right.dimensions, -1)};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_quantity<T2>,
    std::enable_if_t<!std::is_same_v<typename T2::Container, double>, bool>>
CommonQuantityType<T1, T2> fmod(T1 const & left, T2 const & right)
{
    if(right.dimensions == Dimensionless || left.dimensions == right.dimensions)
    {
        return {xt::fmod(left.magnitude, right.magnitude), left.dimensions};
    }
    else
    {
        throw std::runtime_error("Modulo requires same dimensions");
    }
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_quantity<T2>,
    std::enable_if_t<std::is_same_v<typename T2::Container, double>, bool>>
CommonQuantityType<T1, T2> fmod(T1 const & left, T2 const & right)
{
    if(right.dimensions == Dimensionless || left.dimensions == right.dimensions)
    {
        return {std::fmod(left.magnitude, right.magnitude), left.dimensions};
    }
    else
    {
        throw std::runtime_error("Modulo requires same dimensions");
    }
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_not_quantity<T2>,
    std::enable_if_t<!std::is_same_v<std::decay_t<T2>, double>, bool>>
CommonQuantityType<T1, QuantityContainer<T2>>
fmod(T1 const & left, T2 const & right)
{
    return {xt::fmod(left.magnitude, right), left.dimensions};
}

template<
    typename T1, typename T2,
    enable_if_quantity<T1>, enable_if_not_quantity<T2>,
    std::enable_if_t<std::is_same_v<std::decay_t<T2>, double>, bool>>
CommonQuantityType<T1, QuantityContainer<T2>>
fmod(T1 const & left, T2 const & right)
{
    return {std::fmod(left.magnitude, right), left.dimensions};
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
