#include <pybind11/pybind11.h>

#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

void wrap_Species(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::units;

    class_<Species>(m, "Species")
        .def(
            init<Real, Real, Real, Real, Real, Real>(),
            arg("R1"), arg("R2"),
            arg("D")=0, arg("R2_prime")=0, arg("delta_omega")=0, arg("w")=1)
        .def(
            init<Quantity, Quantity, Quantity, Quantity, Quantity, Real>(),
            arg("R1"), arg("R2"),
            arg("D")=0*units::m*units::m/s, arg("R2_prime")=0_Hz,
            arg("delta_omega")=0*rad/s, arg("w")=1)
        .def_readwrite("R1", &Species::R1)
        .def_readwrite("R2", &Species::R2)
        .def_readwrite("D", &Species::D)
        .def_readwrite("R2_prime", &Species::R2_prime)
        .def_readwrite("delta_omega", &Species::delta_omega)
        .def_readwrite("w", &Species::w);
}
