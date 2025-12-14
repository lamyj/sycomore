#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensor.h"
#include "Quantity.h"

void wrap_TensorQ4(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    using TensorQ4 = TensorQ<4>;
    auto TensorQ4Class = wrappers::wrap_quantity_class<TensorQ4>(m, "TensorQ4");
    wrappers::wrap_quantity_array(m, TensorQ4Class);
}
