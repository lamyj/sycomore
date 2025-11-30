#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensorFixed.h"
#include "Quantity.h"

void wrap_Vector4Q(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    auto Vector4QClass = wrappers::wrap_quantity_class<Vector4Q>(m, "Vector4Q");
    wrappers::wrap_quantity_array(Vector4QClass);
}
