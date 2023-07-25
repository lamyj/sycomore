#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include "sycomore/epg/Regular.h"
#include "sycomore/Species.h"

#include "../type_casters.h"

void wrap_epg_Regular(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;

    class_<Regular, Base>(
            m, "Regular",
            "Regular EPG model, where the gradient dephasing is assumed to be "
            "identical during each time interval."
            "\n"
            "In this model, the orders of the model are consecutive positive "
            "integers starting at 0.")
        .def(
            init<
                Species const &, Vector3R const &, unsigned int,
                Quantity const &, double>(),
            "species"_a, "initial_magnetization"_a=Vector3R{0,0,1},
            "initial_size"_a=100, "unit_dephasing"_a=0*units::rad/units::m,
            "gradient_tolerance"_a=1e-5)
        .def(
            init<
                Species const &, Species const &,
                Vector3R const &, Vector3R const &,
                Quantity const &, Quantity const &,
                unsigned int, Quantity const &, double>(),
            "species_a"_a, "species_b"_a, "M0_a"_a, "M0_b"_a, "k_a"_a,
            "delta_b"_a=0*units::Hz, "initial_size"_a=100, 
            "unit_dephasing"_a=0*units::rad/units::m,
            "gradient_tolerance"_a=1e-5)
        .def(
            init<
                Species const &, Quantity const &,
                Vector3R const &, Vector3R const &,
                Quantity const &, 
                unsigned int, Quantity const &, double>(),
            "species_a"_a, "R1_b_or_T1_b"_a, "M0_a"_a, "M0_b"_a, "k_a"_a, 
            "initial_size"_a=100, "unit_dephasing"_a=0*units::rad/units::m,
            "gradient_tolerance"_a=1e-5)
        .def_readwrite("velocity", &Regular::velocity, "Bulk velocity")
        .def_property_readonly(
            "unit_dephasing", &Regular::unit_dephasing,
            "Unit gradient dephasing of the model.")
        .def_property_readonly(
            "orders", &Regular::orders, 
            "The sequence of orders currently stored by the model, in the same "
            "order as the states member. This attribute is read-only.")
        .def(
            "state", overload_cast<std::size_t>(&Base::state, const_),
            "bin"_a,
            "Magnetization at a given state, expressed by its *index*.")
        .def(
            "state", overload_cast<Quantity const &>(&Regular::state, const_),
            "order"_a,
            "Magnetization at a given state, expressed by its *order*.")
        .def(
            "apply_time_interval", 
            static_cast<void(Regular::*)(Quantity const &, Quantity const &)>(
                &Regular::apply_time_interval),
            "duration"_a, "gradient"_a=0*units::T/units::m,
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects.")
        .def(
            "apply_time_interval", 
            static_cast<void(Regular::*)(TimeInterval const &)>(
                &Regular::apply_time_interval),
            "time_interval"_a,
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects.")
        .def(
            "shift", static_cast<void (Regular::*)()>(&Regular::shift), 
            "Apply a unit gradient; in regular EPG, this shifts all orders by 1.")
        .def(
            "shift", 
            static_cast<
                    void (Regular::*)(Quantity const &, Quantity const &)
                >(&Regular::shift), 
            "duration"_a, "gradient"_a,
            "Apply an arbitrary gradient; in regular EPG, this shifts all "
            "orders by an integer number corresponding to a multiple of the "
            "unit gradient.")
        .def(
            "diffusion", &Regular::diffusion, "duration"_a, "gradient"_a,
            "Simulate diffusion during given duration with given gradient "
            "amplitude.")
    ;
}
