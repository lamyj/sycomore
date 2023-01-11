#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include "sycomore/epg/Discrete3D.h"
#include "sycomore/Species.h"

#include "../type_casters.h"

void wrap_epg_Discrete3D(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;

    class_<Discrete3D, Base>(
            m, "Discrete3D",
            "Discrete EPG in which the gradients may be specified in three "
            "dimensions."
        )
        .def(
            init<Species, Vector3R, Quantity>(),
            arg("species"), arg("initial_magnetization")=Vector3R{0,0,1},
            arg("bin_width")=1*units::rad/units::m)
        .def(
            init<
                Species const &, Species const &,
                Vector3R const &, Vector3R const &,
                Quantity const &, Quantity const &, Quantity const &>(),
            arg("species_a"), arg("species_b"), arg("M0_a"), arg("M0_b"),
            arg("k_a"), arg("delta_b")=0*units::Hz,
            arg("bin_width")=1*units::rad/units::m)
        .def(
            init<
                Species const &, Quantity const &,
                Vector3R const &, Vector3R const &,
                Quantity const &, Quantity const &>(),
            arg("species_a"), arg("R1_b_or_T1_b"), arg("M0_a"), arg("M0_b"),
            arg("k_a"), arg("bin_width")=1*units::rad/units::m)
        .def_property_readonly(
            "orders", &Discrete3D::orders, "Orders of the model.")
        .def_property_readonly("bin_width", &Discrete3D::bin_width)
        .def(
            "state", overload_cast<std::size_t>(&Base::state, const_),
            arg("bin"),
            "Magnetization at a given state, expressed by its *index*")
        .def(
            "state",
            overload_cast<Discrete3D::Order const &>(
                &Discrete3D::state, const_),
            arg("order"), "Access a given state of the model")
        .def_property_readonly("elapsed", &Discrete3D::elapsed)
        .def(
            "apply_time_interval",
            overload_cast<Quantity const &, Vector3Q const &>(
                &Discrete3D::apply_time_interval),
            arg("duration"), arg("gradient")=Vector3Q{
                0*units::T/units::m, 0*units::T/units::m, 0*units::T/units::m},
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "apply_time_interval", 
            overload_cast<TimeInterval const &>(
                &Discrete3D::apply_time_interval),
            arg("interval"),
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "shift", &Discrete3D::shift,
            arg("duration"), arg("gradient"),
            "Apply a gradient; in discrete EPG, this shifts all orders by "
            "specified value.")
        .def(
            "diffusion", &Discrete3D::diffusion,
            arg("duration"), arg("gradient"),
            "Simulate diffusion during given duration with given gradient "
            "amplitude.")
    ;
}
