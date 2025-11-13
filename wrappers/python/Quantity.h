#ifndef _d81aa20b_dd7e_4a5a_a47c_837fd494e686
#define _d81aa20b_dd7e_4a5a_a47c_837fd494e686

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "sycomore/Quantity.h"

namespace sycomore
{

namespace wrappers
{

template<typename T>
T as_quantity(pybind11::array_t<pybind11::object> array)
{
    std::vector<size_t> const shape{array.shape(), array.shape()+array.ndim()};
    
    T destination(T::Container::from_shape(shape));
    auto dest_it = destination.magnitude.begin();
    
    array.resize({array.size()});
    
    if(array.size() != 0)
    {
        destination.dimensions = array.data()->cast<sycomore::Quantity>().dimensions;
    }
    
    for(auto && source: array)
    {
        auto const & q = source.cast<sycomore::Quantity>();
        destination.check_dimensions(q, "Constructor requires same dimensions");
        *dest_it = q.magnitude;
        ++dest_it;
    }
    
    array.resize(shape);
    
    return destination;
}

}

}

#endif // _d81aa20b_dd7e_4a5a_a47c_837fd494e686
