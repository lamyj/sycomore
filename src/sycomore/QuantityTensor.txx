#ifndef _962211ec_3f51_400b_923b_6b365e59200a
#define _962211ec_3f51_400b_923b_6b365e59200a

#include "QuantityTensor.h"

#include "sycomore/Quantity.h"

namespace sycomore
{

template<std::size_t N>
TensorQ<N>
::TensorQ(xt::nested_initializer_list_t<double, N> t)
: Base(t)
{
    // Nothing else
}

template<std::size_t N>
TensorQ<N>
::TensorQ(xt::nested_initializer_list_t<Quantity, N> args)
{
    this->_from_array(args);
}

template<std::size_t N>
void
TensorQ<N>
::_from_array(xt::nested_initializer_list_t<Quantity, N> const & args)
{
    auto const dimensions = details::get_dimensions<N>(args);
    details::nested_check_dimensions<N>(args, dimensions);
    
    this->magnitude.resize(xt::shape<shape_type>(args));
    constexpr auto const tmp = xt::layout_type::row_major;
    this->magnitude.layout() == tmp
        ? details::nested_copy_magnitude(this->magnitude.begin(), args)
        : details::nested_copy_magnitude(this->magnitude.template begin<tmp>(), args);
    
    this->dimensions = dimensions;
}

}

namespace std
{

template<std::size_t N>
std::size_t
hash<sycomore::TensorQ<N>>
::operator()(sycomore::TensorQ<N> const & q) const noexcept
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

#endif // _962211ec_3f51_400b_923b_6b365e59200a

