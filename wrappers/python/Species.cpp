#include <algorithm>

#include <pybind11/pybind11.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace
{

void set_D(sycomore::Species & species, pybind11::object const & value)
{
    if(pybind11::isinstance<pybind11::sequence>(value))
    {
        sycomore::Array<sycomore::Quantity> array(pybind11::len(value));
        std::transform(
            value.begin(), value.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        species.set_D(array);
    }
    else
    {
        species.set_D(value.cast<sycomore::Quantity>());
    }
}

sycomore::Species constructor(
    sycomore::Quantity const & R1, sycomore::Quantity const & R2,
    pybind11::object D, sycomore::Quantity const & R2_prime,
    sycomore::Quantity const & delta_omega, sycomore::Real w)
{
    sycomore::Species species(
        R1, R2, {0, sycomore::Diffusion}, R2_prime, delta_omega, w);
    set_D(species, D);
    return species;
}

}

void wrap_Species(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::units;

    class_<Species>(m, "Species")
        .def(
            init(&constructor),
            arg("R1"), arg("R2"),
            arg("D")=0*units::m*units::m/s, arg("R2_prime")=0_Hz,
            arg("delta_omega")=0*rad/s, arg("w")=1)
        .def_property("R1", &Species::get_R1, &Species::set_R1)
        .def_property_readonly("T1", &Species::get_T1)
        .def_property("R2", &Species::get_R2, &Species::set_R2)
        .def_property_readonly("T2", &Species::get_T2)
        .def_property("D", &Species::get_D, set_D)
        .def_property(
            "R2_prime", &Species::get_R2_prime, &Species::set_R2_prime)
        .def_property_readonly("T2_prime", &Species::get_T2_prime)
        .def_property(
            "delta_omega", &Species::get_delta_omega, &Species::set_delta_omega)
        .def_readwrite("w", &Species::w);
}
