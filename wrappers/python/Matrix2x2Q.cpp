#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensorFixed.h"
#include "Quantity.h"

void wrap_Matrix2x2Q(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    auto Matrix2x2QClass = wrappers::wrap_quantity_class<Matrix2x2Q>(m, "Matrix2x2Q");
    wrappers::wrap_quantity_array(m, Matrix2x2QClass);
}
