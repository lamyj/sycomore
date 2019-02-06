#include <pybind11/pybind11.h>

void wrap_Dimensions(pybind11::module &);
void wrap_Quantity(pybind11::module &);
void wrap_units(pybind11::module &);

void wrap_magnetization(pybind11::module &);

void wrap_GridScanner(pybind11::module &);
void wrap_Array(pybind11::module &);
void wrap_Grid(pybind11::module &);

void wrap_Pulse(pybind11::module &);
void wrap_Species(pybind11::module &);
void wrap_TimeInterval(pybind11::module &);

void wrap_Model(pybind11::module &);

PYBIND11_MODULE(_sycomore, _sycomore)
{
    wrap_Dimensions(_sycomore);
    wrap_Quantity(_sycomore);
    wrap_units(_sycomore);

    wrap_magnetization(_sycomore);

    wrap_GridScanner(_sycomore);
    wrap_Array(_sycomore);
    wrap_Grid(_sycomore);

    wrap_Pulse(_sycomore);
    wrap_Species(_sycomore);
    wrap_TimeInterval(_sycomore);

    wrap_Model(_sycomore);
}
