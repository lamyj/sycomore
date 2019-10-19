#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "sycomore/epg/Discrete3D.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"

void wrap_epg_Discrete3D(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;

    class_<Discrete3D>(m, "Discrete3D")
        .def(
            init<Species, Magnetization, Quantity>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("bin_width")=1*units::rad/units::m)
        .def_readwrite("species", &Discrete3D::species)
        .def_property_readonly("orders", &Discrete3D::orders)
        .def_property_readonly("bin_width", &Discrete3D::bin_width)
        .def("state", &Discrete3D::state, arg("order"))
        .def_property_readonly(
            "states", [](Discrete3D const & model){
                auto const states_cpp = model.states();
                std::vector<std::size_t> const shape{model.size(),3};
                auto && data = states_cpp.data();
                array_t<Complex> states_py(shape, data);
                return states_py;
            })
        .def_property_readonly("echo", &Discrete3D::echo)
        .def(
            "apply_pulse", &Discrete3D::apply_pulse,
            arg("angle"), arg("phase")=0*units::rad)
        .def(
            "apply_time_interval",
            [](
                Discrete3D & model, Quantity const & duration,
                sequence const & gradient, Real threshold)
            {
                Array<Quantity> array(gradient.size());
                for(std::size_t i=0; i<array.size(); ++i)
                {
                    array[i] = gradient[i].cast<Quantity>();
                }
                model.apply_time_interval(duration, array, threshold);
            },
            arg("duration"), arg("gradient")=Array<Quantity>{
                0*units::T/units::m, 0*units::T/units::m, 0*units::T/units::m},
            arg("threshold")=0.)
        .def("shift", &Discrete3D::shift, arg("duration"), arg("gradient"))
        .def("relaxation", &Discrete3D::relaxation, arg("duration"))
        .def("diffusion", &Discrete3D::diffusion, arg("duration"), arg("gradient"))
        .def("__len__", &Discrete3D::size)
    ;
}
