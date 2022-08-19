#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

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
            init<Species, Magnetization, unsigned int, Quantity const &, double>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("initial_size")=100, 
            arg("unit_gradient_area")=0*units::mT/units::m*units::ms, 
            arg("gradient_tolerance")=1e-5)
        .def_property(
            "species",
            [](Regular const & r){ return r.get_species();},
            [](Regular & r, Species const & s){ return r.set_species(s);})
        .def_readwrite("threshold", &Regular::threshold)
        .def_readwrite("delta_omega", &Regular::delta_omega)
        .def_readwrite("velocity", &Regular::velocity)
        .def_property_readonly(
            "unit_gradient_area", &Regular::unit_gradient_area,
            "Unit gradient area of the model.")
        .def_property_readonly(
            "orders", &Regular::orders, 
            "The sequence of orders currently stored by the model, in the same "
            "order as the states member. This attribute is read-only.")
        .def(
            "state",
            static_cast<
                    std::vector<Complex> (Regular::*)(std::size_t) const
                >(&Regular::state),
            arg("bin"),
            "Magnetization at a given state, expressed by its *index*")
        .def(
            "state",
            static_cast<
                    std::vector<Complex> (Regular::*)(Quantity const &) const
                >(&Regular::state),
            arg("order"),
            "Magnetization at a given state, expressed by its *order*.")
        .def_property_readonly(
            "states", [](Regular const & model){
                auto const states_cpp = model.states();
                std::vector<std::size_t> const shape{model.size(), 3};
                Complex const * data = states_cpp.data();
                array_t<Complex> states_py(shape, data);
                return states_py;
            },
            "Return all states in the model, where each state is stored as "
            "F_k, F*_{-k}, Z_k, in order of increasing order.")
        .def_property_readonly(
            "echo", [](Regular const & r){ return r.echo(); },
            "Echo signal, i.e. F_0")
        .def(
            "apply_pulse", 
            static_cast<void(Regular::*)(Quantity const &, Quantity const &)>(
                &Regular::apply_pulse),
            arg("angle"), arg("phase")=0*units::rad,
            "Apply an RF hard pulse.")
        .def(
            "apply_time_interval", 
            static_cast<void(Regular::*)(Quantity const &, Quantity const &)>(
                &Regular::apply_time_interval),
            arg("duration"), arg("gradient")=0*units::T/units::m,
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects.")
        .def(
            "apply_time_interval", 
            static_cast<void(Regular::*)(TimeInterval const &)>(
                &Regular::apply_time_interval),
            arg("time_interval"),
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
            arg("duration"), arg("gradient"),
            "Apply an arbitrary gradient; in regular EPG, this shifts all "
            "orders by an integer number corresponding to a multiple of the "
            "unit gradient.")
        .def(
            "relaxation", &Regular::relaxation, arg("duration"),
            "Simulate the relaxation during given duration.")
        .def(
            "diffusion", &Regular::diffusion, arg("duration"), arg("gradient"),
            "Simulate diffusion during given duration with given gradient "
            "amplitude.")
        .def(
            "off_resonance", &Regular::off_resonance, 
            arg("duration"),
            "Simulate field- and species related off-resonance effects during "
            "given duration with given frequency offset.")
        .def("__len__", &Regular::size, "Number of states of the model")
    ;
}
