#include <pybind11/pybind11.h>

#include "sycomore/Pulse.h"
#include "sycomore/Quantity.h"

#include "type_casters.h"

void wrap_Pulse(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<Pulse>(m, "Pulse", "RF pulse")
        .def(init<Quantity, Quantity>(), arg("angle"), arg("phase")=0*units::rad)
        .def_property(
            "angle", &Pulse::angle, &Pulse::set_angle,
            "Flip angle of the pulse")
        .def_property(
            "phase", &Pulse::phase, &Pulse::set_phase, "Phase of the pulse");
}
