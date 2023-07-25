#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include "sycomore/epg/Discrete.h"
#include "sycomore/Species.h"

#include "../type_casters.h"

void wrap_epg_Discrete(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;
    
    class_<Discrete, Base>(m, "Discrete", 
        "Discrete EPG model, where the gradient dephasing may vary across time "
        "intervals."
        "\n"
        "In this model, the orders of the model are stored in bins of "
        "user-specified width (hence the term \"discrete\"), expressed in "
        "rad/m.")
        .def(
            init<Species const &, Vector3R const &, Quantity const &>(),
            "species"_a, "initial_magnetization"_a=Vector3R{0,0,1},
            "bin_width"_a=1*units::rad/units::m)
        .def(
            init<
                Species const &, Species const &,
                Vector3R const &, Vector3R const &,
                Quantity const &, Quantity const &, Quantity const &>(),
            "species_a"_a, "species_b"_a, "M0_a"_a, "M0_b"_a, "k_a"_a,
            "delta_b"_a=0*units::Hz, "bin_width"_a=1*units::rad/units::m)
        .def(
            init<
                Species const &, Quantity const &,
                Vector3R const &, Vector3R const &,
                Quantity const &, Quantity const &>(),
            "species_a"_a, "R1_b_or_T1_b"_a, "M0_a"_a, "M0_b"_a, "k_a"_a,
            "bin_width"_a=1*units::rad/units::m)
        .def_property_readonly(
            "orders", &Discrete::orders, 
            "The sequence of orders currently stored by the model, in the same "
            "order as the states member. This attribute is read-only.")
        .def_property_readonly("bin_width", &Discrete::bin_width)
        .def(
            "state", overload_cast<std::size_t>(&Base::state, const_),
            "bin"_a,
            "Magnetization at a given state, expressed by its *index*")
        .def(
            "state", overload_cast<Quantity const &>(&Discrete::state, const_),
            "order"_a,
            "Magnetization at a given state, expressed by its *order*.")
        .def(
            "apply_time_interval", 
            static_cast<void(Discrete::*)(Quantity const &, Quantity const &)>(
                &Discrete::apply_time_interval),
            "duration"_a, "gradient"_a=0*units::T/units::m,
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "apply_time_interval", 
            static_cast<void(Discrete::*)(TimeInterval const &)>(
                &Discrete::apply_time_interval),
            "interval"_a,
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "shift", &Discrete::shift, "duration"_a, "gradient"_a,
            "Apply a gradient; in discrete EPG, this shifts all orders by" 
            "specified value.")
        .def(
            "diffusion", &Discrete::diffusion, "duration"_a, "gradient"_a,
            "Simulate diffusion during given duration with given gradient ",
            "amplitude.")
    ;
}
