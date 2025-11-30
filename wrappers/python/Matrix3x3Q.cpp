#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensorFixed.h"
#include "Quantity.h"

void wrap_Matrix3x3Q(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    auto Matrix3x3QClass = wrappers::wrap_quantity_class<Matrix3x3Q>(m, "Matrix3x3Q");
    wrappers::wrap_quantity_array(Matrix3x3QClass);
}
