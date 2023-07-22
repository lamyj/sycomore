#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "sycomore/epg/Base.h"
#include "sycomore/Species.h"

#include "../type_casters.h"

void wrap_epg_Base(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    using namespace sycomore;
    using namespace sycomore::epg;

    class_<Base>(m, "Base", "Base class for all EPG models")
        // Pure virtual "size" function -> no constructor
        .def_property(
            "species",
            [](Base const & b){ return b.species();},
            [](Base & b, Species const & s){ return b.set_species(s);})
        .def("M0", &Base::M0, "pool"_a=0)
        .def(
            "set_M0", overload_cast<std::size_t, Real const &>(&Base::set_M0),
            "pool"_a, "M0"_a)
        .def("set_M0", overload_cast<Real const &>(&Base::set_M0), "M0"_a)
        .def("k", &Base::k, "pool"_a)
        .def("set_k", &Base::set_k, "pool"_a, "M0"_a)
        .def_property(
            "delta_b",
            [](Base const & b){ return b.delta_b();},
            [](Base & b, Quantity const & q){ return b.set_delta_b(q);})
        .def_readwrite("threshold", &Base::threshold)
        .def_readwrite("delta_omega", &Base::delta_omega)
        .def_property_readonly("kind", &Base::kind)
        .def_property_readonly("pools", &Base::pools)
        .def_property_readonly(
            "states", &Base::states,
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
