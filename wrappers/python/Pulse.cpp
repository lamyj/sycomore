#include <pybind11/pybind11.h>

#include "sycomore/Pulse.h"
#include "sycomore/Quantity.h"

void wrap_Pulse(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<Pulse>(m, "Pulse")
        .def(init<Quantity, Quantity>())
        .def_property("angle", &Pulse::get_angle, &Pulse::set_angle)
        .def_property("phase", &Pulse::get_phase, &Pulse::set_phase);
}
