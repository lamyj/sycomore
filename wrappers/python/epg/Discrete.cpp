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
    
    class_<Discrete>(m, "Discrete", 
        "Discrete EPG model, where the gradient moments may vary across time "
        "intervals."
        "\n"
        "In this model, the orders of the model are stored in bins of "
        "user-specified width (hence the term \"discrete\"), expressed in "
        "rad/m.")
        .def(
            init<Species, Magnetization, Quantity>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("bin_width")=1*units::rad/units::m)
        .def_readwrite("species", &Discrete::species)
        .def_readwrite("delta_omega", &Discrete::delta_omega)
        .def_readwrite("threshold", &Discrete::threshold)
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
        .def_property_readonly(
            "states", [](Discrete const & model){
                auto const states_cpp = model.states();
                std::vector<std::size_t> const shape{model.orders().size(),3};
                auto && data = states_cpp.data();
                array_t<Complex> states_py(shape, data);
                return states_py;
            },
            "The sequence of states currently stored by the model, in the same "
            "order as the orders member. This attribute is a read-only, 3×N "
            "array of complex numbers.")
        .def_property_readonly(
            "echo", &Discrete::echo, "The echo signal, i.e. F_0.")
        .def(
            "apply_pulse", &Discrete::apply_pulse,
            arg("angle"), arg("phase")=0*units::rad,
            "Apply an RF hard pulse.")
        .def(
            "apply_time_interval", 
            static_cast<void(Discrete::*)(Quantity const &, Quantity const &, Real)>(
                &Discrete::apply_time_interval),
            arg("duration"), arg("gradient")=0*units::T/units::m,
            arg("threshold")=0.,
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
            "relaxation", &Discrete::relaxation, arg("duration"),
            "Simulate the relaxation during given duration.")
        .def(
            "diffusion", &Discrete::diffusion, arg("duration"), arg("gradient"),
            "Simulate diffusion during given duration with given gradient ",
            "amplitude.")
        .def(
            "off_resonance", &Discrete::off_resonance, 
            arg("duration"),
            "Simulate field- and species related off-resonance effects during "
            "given duration with given frequency offset.")
        .def("__len__", &Discrete::size, "Number of states of the model.")
    ;
}
