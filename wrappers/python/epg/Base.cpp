#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "sycomore/epg/Base.h"
#include "sycomore/Species.h"

void wrap_epg_Base(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;

    class_<Base>(m, "Base", "Base class for all EPG models")
        // Pure virtual "size" function -> no constructor
        .def_property(
            "species",
            [](Base const & r){ return r.get_species();},
            [](Base & r, Species const & s){ return r.set_species(s);})
        .def_readwrite("threshold", &Base::threshold)
        .def_readwrite("delta_omega", &Base::delta_omega)
        .def_property_readonly("kind", &Base::kind)
        .def_property_readonly("pools", &Base::pools)
        .def_property_readonly(
            "states", [](Base const & model){
                std::vector<std::size_t> shape;
                if(model.kind() == Model::SinglePool)
                {
                    shape = {model.size(), 3};
                }
                else 
                {
                    shape = {model.size(), model.pools(), 3};
                }
                
                auto const states_cpp = model.states();
                auto const * data = states_cpp.data();
                
                array_t<Complex> states_py(shape, data);
                return states_py;
            },
            "Return all states in the model, where each state is stored as "
            "F_k, F*_{-k}, Z_k, in order of increasing order.")
        .def_property_readonly("elapsed", &Base::elapsed)
        .def_property_readonly(
            "echo", [](Base const & r){ return r.echo(); },
            "Echo signal, i.e. F_0")
        .def(
            "apply_pulse", 
            static_cast<void(Base::*)(Quantity const &, Quantity const &)>(
                &Base::apply_pulse),
            arg("angle"), arg("phase")=0*units::rad,
            "Apply an RF hard pulse.")
        .def(
            "apply_pulse", 
            static_cast<
                    void(Base::*)(
                        Quantity const &, Quantity const &,
                        Quantity const &, Quantity const &)
                >(&Base::apply_pulse),
            arg("angle_a"), arg("phase_a"), arg("angle_B"), arg("phase_b"),
            "Apply an RF hard pulse.")
        .def(
            "apply_pulse", 
            static_cast<
                    void(Base::*)(Quantity const &, Quantity const &, Real)
                >(&Base::apply_pulse),
            arg("angle_a"), arg("phase_a"), arg("saturation"),
            "Apply an RF hard pulse.")
        .def(
            "relaxation", &Base::relaxation, arg("duration"),
            "Simulate the relaxation during given duration.")
        .def(
            "off_resonance", &Base::off_resonance, 
            arg("duration"),
            "Simulate field- and species related off-resonance effects during "
            "given duration with given frequency offset.")
        .def("__len__", &Base::size, "Number of states of the model")
    ;
}
