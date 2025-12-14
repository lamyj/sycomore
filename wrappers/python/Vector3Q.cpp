#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensorFixed.h"
#include "Quantity.h"

void wrap_Vector3Q(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    auto Vector3QClass = wrappers::wrap_quantity_class<Vector3Q>(m, "Vector3Q");
    wrappers::wrap_quantity_array(m, Vector3QClass);
}
