#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

void wrap_Dimensions(pybind11::module &);
void wrap_Quantity(pybind11::module &);
void wrap_units(pybind11::module &);

PYBIND11_MODULE(_sycomore, _sycomore)
{
    wrap_Dimensions(_sycomore);
    wrap_Quantity(_sycomore);
    wrap_units(_sycomore);
}
