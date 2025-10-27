#ifndef _c268e431_674b_4c07_a919_3f7c4d4f3273
#define _c268e431_674b_4c07_a919_3f7c4d4f3273

#include "QuantityTensorFixed.h"

#include "sycomore/Quantity.h"

namespace sycomore
{

template<typename S>
TensorFixedQ<S>
::TensorFixedQ(xt::nested_initializer_list_t<double, std::tuple_size<S>::value> t)
: Base(t)
{
    // Nothing else
}

template<typename S>
TensorFixedQ<S>
::TensorFixedQ(xt::nested_initializer_list_t<Quantity, rank> args)
{
    this->_from_array(args);
}

template<typename S>
void
TensorFixedQ<S>
::_from_array(xt::nested_initializer_list_t<Quantity, rank> const & args)
{
    auto const dimensions = details::get_dimensions<rank>(args);
    details::nested_check_dimensions<rank>(args, dimensions);
    
    constexpr auto const tmp = xt::layout_type::row_major;
    this->magnitude.layout() == tmp
        ? details::nested_copy_magnitude(this->magnitude.begin(), args)
        : details::nested_copy_magnitude(this->magnitude.template begin<tmp>(), args);
    
    this->dimensions = dimensions;
}

}

namespace std
{
    
template<typename S>
std::size_t
hash<sycomore::TensorFixedQ<S>>
::operator()(sycomore::TensorFixedQ<S> const & q) const noexcept
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

#endif // _c268e431_674b_4c07_a919_3f7c4d4f3273
