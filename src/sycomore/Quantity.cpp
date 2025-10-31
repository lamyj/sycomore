#include "sycomore/Quantity.h"

#include "sycomore/QuantityBase.h"

namespace sycomore
{

Quantity
::Quantity(QuantityConstReference const & q)
: Base(q.magnitude, q.dimensions)
{
    // Nothing else
}

template<template<typename> typename Operator>
bool order(Quantity const & left, Quantity const & right)
{
    left.check_dimensions(right, "Ordering requires same dimensions");
    return Operator<Quantity::Container>()(left.magnitude, right.magnitude);
}

bool operator<(Quantity const & left, Quantity const & right)
{
    return order<std::less>(left, right);
}

bool operator<=(Quantity const & left, Quantity const & right)
{
    return order<std::less_equal>(left, right);
}

bool operator>(Quantity const & left, Quantity const & right)
{
    return order<std::greater>(left, right);
}

bool operator>=(Quantity const & left, Quantity const & right)
{
    return order<std::greater_equal>(left, right);
}

namespace details
{

template<>
Dimensions get_dimensions<1>(xt::nested_initializer_list_t<Quantity, 1> const & t)
{
    return t.begin()->dimensions;
}

template<>
void nested_check_dimensions<1>(
    xt::nested_initializer_list_t<Quantity, 1> const & t,
    Dimensions const & dimensions)
{
    for(auto && x: t)
    {
        if(x.dimensions != dimensions)
        {
            throw std::runtime_error("Constructor requires same dimensions");
        }
    }
}

}

}

namespace std
{

std::size_t
hash<sycomore::Quantity>
::operator()(sycomore::Quantity const & q) const noexcept
{
    std::size_t seed=0;
    hash<double> hasher;
    sycomore::combine_hashes(seed, hasher(q.magnitude));
    sycomore::combine_hashes(
        seed, std::hash<sycomore::Dimensions>{}(q.dimensions));
    return seed;
}

}
