#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensorFixed.h"
#include "Quantity.h"

void wrap_Vector2Q(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    auto Vector2QClass = wrappers::wrap_quantity_class<Vector2Q>(m, "Vector2Q");
    wrappers::wrap_quantity_array(Vector2QClass);
}
