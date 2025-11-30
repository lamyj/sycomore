#include <pybind11/pybind11.h>

#include "sycomore/QuantityTensor.h"
#include "Quantity.h"

void wrap_TensorQ1(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    
    using TensorQ1 = TensorQ<1>;
    auto TensorQ1Class = wrappers::wrap_quantity_class<TensorQ1>(m, "TensorQ1");
    wrappers::wrap_quantity_array(TensorQ1Class);
}
