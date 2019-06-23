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

    class_<Regular>(m, "Regular")
        .def(
            init<Species, Magnetization, unsigned int>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("initial_size")=100)
        .def_readwrite("species", &Regular::species)
        .def_property_readonly("states_count", &Regular::states_count)
        .def("state", &Regular::state, arg("order"))
        .def_property_readonly(
            "states", [](Regular const & model){
                auto const states_cpp = model.states();
                std::vector<std::size_t> const shape{model.states_count(),3};
                Complex const * data = states_cpp.data();
                array_t<Complex> states_py(shape, data);
                return states_py;
            })
        .def_property_readonly("echo", &Regular::echo)
        .def(
            "apply_pulse", &Regular::apply_pulse,
            arg("angle"), arg("phase")=0*units::rad)
        .def(
            "apply_time_interval", &Regular::apply_time_interval,
            arg("duration"), arg("gradient")=0*units::T/units::m)
        .def("shift", &Regular::shift)
        .def("relaxation", &Regular::relaxation, arg("duration"))
        .def("diffusion", &Regular::diffusion, arg("duration"), arg("gradient"))
    ;
}
