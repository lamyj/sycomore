#ifndef _3be9ab71_0f7d_4be8_84d2_2c6b925d98e1
#define _3be9ab71_0f7d_4be8_84d2_2c6b925d98e1

#include "QuantityArray.h"

#include "sycomore/Quantity.h"

namespace sycomore
{

template<std::size_t D>
void
ArrayQ
::_from_array(xt::nested_initializer_list_t<Quantity, D> const & args)
{
    auto const dimensions = details::get_dimensions<D>(args);
    details::nested_check_dimensions<D>(args, dimensions);
    
    this->magnitude.resize(xt::shape<shape_type>(args));
    constexpr auto const tmp = xt::layout_type::row_major;
    this->magnitude.layout() == tmp
        ? details::nested_copy_magnitude(this->magnitude.begin(), args)
        : details::nested_copy_magnitude(this->magnitude.begin<tmp>(), args);
    
    this->dimensions = dimensions;
}

}

#endif // _3be9ab71_0f7d_4be8_84d2_2c6b925d98e1
