#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "sycomore/epg/Discrete.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"

void wrap_epg_Discrete(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;
    
    class_<Discrete, Base>(m, "Discrete", 
        "Discrete EPG model, where the gradient moments may vary across time "
        "intervals."
        "\n"
        "In this model, the orders of the model are stored in bins of "
        "user-specified width (hence the term \"discrete\"), expressed in "
        "rad/m.")
        .def(
            init<Species const &, Magnetization const &, Quantity const &>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("bin_width")=1*units::rad/units::m)
        .def(
            init<
                Species const &, Species const &,
                Magnetization const &, Magnetization const &,
                Quantity const &, Quantity const &, Quantity const &>(),
            arg("species_a"), arg("species_b"), arg("M0_a"), arg("M0_b"),
            arg("k_a"), arg("delta_b")=0*units::Hz,
            arg("bin_width")=1*units::rad/units::m)
        .def(
            init<
                Species const &, Quantity const &,
                Magnetization const &, Magnetization const &,
                Quantity const &, Quantity const &>(),
            arg("species_a"), arg("R1_b_or_T1_b"), arg("M0_a"), arg("M0_b"),
            arg("k_a"), arg("bin_width")=1*units::rad/units::m)
        .def_property_readonly(
            "orders", &Discrete::orders, 
            "The sequence of orders currently stored by the model, in the same "
            "order as the states member. This attribute is read-only.")
        .def_property_readonly("bin_width", &Discrete::bin_width)
        .def(
            "state",
            static_cast<
                    std::vector<Complex> (Discrete::*)(std::size_t) const
                >(&Discrete::state),
            arg("bin"),
            "Magnetization at a given state, expressed by its *index*")
        .def(
            "state",
            static_cast<
                    std::vector<Complex> (Discrete::*)(Quantity const &) const
                >(&Discrete::state),
            arg("order"),
            "Magnetization at a given state, expressed by its *order*.")
        .def(
            "apply_time_interval", 
            static_cast<void(Discrete::*)(Quantity const &, Quantity const &)>(
                &Discrete::apply_time_interval),
            arg("duration"), arg("gradient")=0*units::T/units::m,
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "apply_time_interval", 
            static_cast<void(Discrete::*)(TimeInterval const &)>(
                &Discrete::apply_time_interval),
            arg("interval"),
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "shift", &Discrete::shift, arg("duration"), arg("gradient"),
            "Apply a gradient; in discrete EPG, this shifts all orders by" 
            "specified value.")
        .def(
            "diffusion", &Discrete::diffusion, arg("duration"), arg("gradient"),
            "Simulate diffusion during given duration with given gradient ",
            "amplitude.")
    ;
}
