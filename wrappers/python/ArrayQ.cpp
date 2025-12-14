#include <pybind11/pybind11.h>

#include "sycomore/QuantityArray.h"
#include "Quantity.h"

void wrap_ArrayQ(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    auto ArrayQClass = wrappers::wrap_quantity_class<ArrayQ>(m, "ArrayQ");
    wrappers::wrap_quantity_array(m, ArrayQClass);
}
