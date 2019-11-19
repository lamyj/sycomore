#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "sycomore/epg/Regular.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"

void wrap_epg_Regular(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;

    class_<Regular>(
            m, "Regular",
            "Regular EPG model, where the gradient moment is assumed to be "
            "identical during each time interval."
            "\n"
            "In this model, the orders of the model are consecutive positive "
            "integers starting at 0.")
        .def(
            init<Species, Magnetization, unsigned int>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("initial_size")=100)
        .def_readwrite("species", &Regular::species)
        .def_property_readonly(
            "states_count", &Regular::states_count, 
            "Number of states in the model.")
        .def(
            "state", &Regular::state, arg("order"), 
            "Access a given state of the model.")
        .def_property_readonly(
            "states", [](Regular const & model){
                auto const states_cpp = model.states();
                std::vector<std::size_t> const shape{model.states_count(),3};
                Complex const * data = states_cpp.data();
                array_t<Complex> states_py(shape, data);
                return states_py;
            },
            "Return all states in the model, where each state is stored as "
            "F̃_k, F̃^*_{-k}, Z̃_k, in order of increasing order.")
        .def_property_readonly("echo", &Regular::echo, "Echo signal, i.e. F̃_0")
        .def(
            "apply_pulse", &Regular::apply_pulse,
            arg("angle"), arg("phase")=0*units::rad,
            "Apply an RF hard pulse.")
        .def(
            "apply_time_interval", &Regular::apply_time_interval,
            arg("duration"), arg("gradient")=0*units::T/units::m,
            "Apply a time interval, i.e. relaxation, diffusion, and gradient.")
        .def(
            "shift", &Regular::shift, 
            "Apply a gradient; in regular EPG, this shifts all orders by 1.")
        .def(
            "relaxation", &Regular::relaxation, arg("duration"),
            "Simulate the relaxation during given duration.")
        .def(
            "diffusion", &Regular::diffusion, arg("duration"), arg("gradient"),
            "Simulate diffusion during given duration with given gradient "
            "amplitude.")
    ;
}
