#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensor.h"
#include "Quantity.h"

void wrap_TensorQ2(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    using TensorQ2 = TensorQ<2>;
    auto TensorQ2Class = wrappers::wrap_quantity_class<TensorQ2>(m, "TensorQ2");
    wrappers::wrap_quantity_array(TensorQ2Class);
}
