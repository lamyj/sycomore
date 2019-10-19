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

    class_<Discrete>(m, "Discrete")
        .def(
            init<Species, Magnetization, Quantity>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("bin_width")=1*units::rad/units::m)
        .def_readwrite("species", &Discrete::species)
        .def_property_readonly("orders", &Discrete::orders)
        .def_property_readonly("bin_width", &Discrete::bin_width)
        .def(
            "state",
            static_cast<
                    std::vector<Complex> (Discrete::*)(std::size_t) const
                >(&Discrete::state),
            arg("bin"))
        .def(
            "state",
            static_cast<
                    std::vector<Complex> (Discrete::*)(Quantity const &) const
                >(&Discrete::state),
            arg("order"))
        .def_property_readonly(
            "states", [](Discrete const & model){
                auto const states_cpp = model.states();
                std::vector<std::size_t> const shape{model.orders().size(),3};
                auto && data = states_cpp.data();
                array_t<Complex> states_py(shape, data);
                return states_py;
            })
        .def_property_readonly("echo", &Discrete::echo)
        .def(
            "apply_pulse", &Discrete::apply_pulse,
            arg("angle"), arg("phase")=0*units::rad)
        .def(
            "apply_time_interval", &Discrete::apply_time_interval,
            arg("duration"), arg("gradient")=0*units::T/units::m,
            arg("threshold")=0.)
        .def("shift", &Discrete::shift, arg("duration"), arg("gradient"))
        .def("relaxation", &Discrete::relaxation, arg("duration"))
        .def("diffusion", &Discrete::diffusion, arg("duration"), arg("gradient"))
        .def("__len__", &Discrete::size)
    ;
}
