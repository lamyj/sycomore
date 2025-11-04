#include "QuantityArray.h"

#include "sycomore/Quantity.h"

namespace sycomore
{

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<double, 1> t)
: Base(t)
{
    // Nothing else
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<Quantity, 1> args)
{
    this->_from_array<1>(args);
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<double, 2> t)
: Base(t)
{
    // Nothing else
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<Quantity, 2> args)
{
    this->_from_array<2>(args);
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<double, 3> t)
: Base(t)
{
    // Nothing else
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<Quantity, 3> args)
{
    this->_from_array<3>(args);
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<double, 4> t)
: Base(t)
{
    // Nothing else
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<Quantity, 4> args)
{
    this->_from_array<4>(args);
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<double, 5> t)
: Base(t)
{
    // Nothing else
}

ArrayQ
::ArrayQ(xt::nested_initializer_list_t<Quantity, 5> args)
{
    this->_from_array<5>(args);
}

}

namespace std
{

std::size_t
hash<sycomore::ArrayQ>
::operator()(sycomore::ArrayQ const & q) const noexcept
{
    std::size_t seed=0;
    hash<double> hasher;
    for(auto && x: q.magnitude)
    {
        sycomore::combine_hashes(seed, hasher(x));
    }
    sycomore::combine_hashes(
        seed, std::hash<sycomore::Dimensions>{}(q.dimensions));
    return seed;
}

}
