#include <pybind11/pybind11.h>

#include "sycomore/Pulse.h"
#include "sycomore/Quantity.h"

void wrap_Pulse(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<Pulse>(m, "Pulse")
        .def(init<Real, Real>())
        .def(init<Quantity, Quantity>())
        .def_readwrite("angle", &Pulse::angle)
        .def_readwrite("phase", &Pulse::phase);
}
