#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensor.h"
#include "Quantity.h"

void wrap_TensorQ3(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    using TensorQ3 = TensorQ<3>;
    auto TensorQ3Class = wrappers::wrap_quantity_class<TensorQ3>(m, "TensorQ3");
    wrappers::wrap_quantity_array(TensorQ3Class);
}
