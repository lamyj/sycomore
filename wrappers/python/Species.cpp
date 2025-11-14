#include <algorithm>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace
{

void set_D(sycomore::Species & species, pybind11::object const & value)
{
    try
    {
        species.set_D(value.cast<sycomore::Quantity>());
    }
    catch(pybind11::cast_error const &)
    {
        species.set_D(value.cast<sycomore::Matrix3x3Q>());
    }
}

}

void wrap_Species(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::units;

    class_<Species>(
            m, "Species", "Species described by its NMR parameters")
        .def(
            init<
                Quantity const &, Quantity const &,
                Quantity const &, Quantity const &>(),
            "R1"_a, "R2"_a, "D"_a=0*units::m*units::m/units::s,
            "delta_omega"_a=0*units::Hz)
        .def(
            init<
                Quantity const &, Quantity const &,
                Matrix3x3Q const &, Quantity const &>(),
            "R1"_a, "R2"_a, "D"_a, "delta_omega"_a=0*units::Hz)
        .def_property(
            "R1", &Species::R1, &Species::set_R1,
            "Longitudinal relaxation rate.")
        .def_property(
            "T1", &Species::T1, &Species::set_T1,
            "Longitudinal relaxation time.")
        .def_property(
            "R2", &Species::R2, &Species::set_R2, 
            "Transversal relaxation rate.")
        .def_property(
            "T2", &Species::T2, &Species::set_T2,
            "Transversal relaxation time.")
        .def_property("D", &Species::D, set_D, "Diffusion tensor.")
        .def_property(
            "delta_omega", &Species::delta_omega, &Species::set_delta_omega,
            "Frequency offset.")
        .def(pickle(
            [](Species const & s) {
                return make_tuple(s.R1(), s.R2(), s.D(), s.delta_omega());
            },
            [](tuple t) {
                if(t.size() != 4)
                {
                    throw std::runtime_error("Invalid state!");
                }
                return Species(
                    t[0].cast<Quantity>(), t[1].cast<Quantity>(),
                    t[2].cast<Matrix3x3Q>(), t[3].cast<Quantity>());
            }
        ));
}
