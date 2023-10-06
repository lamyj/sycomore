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
        .def(
            "M0", &Base::M0, "pool"_a=0,
            "Return the equilibrium magnetization of one of the pools")
        .def(
            "set_M0", overload_cast<std::size_t, Real const &>(&Base::set_M0),
            "pool"_a, "M0"_a,
            "Set the equilibrium magnetization of one of the pools")
        .def(
            "set_M0", overload_cast<Real const &>(&Base::set_M0), "M0"_a,
            "Set the equilibrium magnetization of one of the pools")
        .def(
            "k", &Base::k, "pool"_a,
            "Return the exchange constant of one of the pools")
        .def(
            "set_k", &Base::set_k, "pool"_a, "M0"_a,
            "Set the exchange constant of one of the pools")
        .def_property(
            "delta_b",
            [](Base const & b){ return b.delta_b();},
            [](Base & b, Quantity const & q){ return b.set_delta_b(q);},
            "Frequency offset of the second pool")
        .def_readwrite(
            "threshold", &Base::threshold,
            "Threshold used to cull states with low population")
        .def_readwrite(
            "delta_omega", &Base::delta_omega,
            "Frequency offset of the simulator")
        .def_property_readonly(
            "kind", &Base::kind,
            "Return the kind of the model, set at creation")
        .def_property_readonly(
            "pools", &Base::pools,
            "Return the number of pools of the model, set at creation")
        .def_property_readonly(
            "states", &Base::states,
            R"(Return all states in the model, as :math:`\tilde{F}`, )"
                R"(:math:`\tilde{F}^*`, and :math:`\tilde{Z}` for each order )"
                "and each pool.")
        .def_property_readonly(
            "elapsed", &Base::elapsed, "Return the elapsed time")
        .def_property_readonly(
            "echo", [](Base const & r){ return r.echo(); },
            "Echo signal, i.e. :math:`F_0`")
        .def(
            "apply_pulse", 
            static_cast<void(Base::*)(Quantity const &, Quantity const &)>(
                &Base::apply_pulse),
            "angle"_a, "phase"_a=0*units::rad,
            "Apply an RF hard pulse to a single-pool model")
        .def(
            "apply_pulse", 
            static_cast<
                    void(Base::*)(
                        Quantity const &, Quantity const &,
                        Quantity const &, Quantity const &)
                >(&Base::apply_pulse),
            "angle_a"_a, "phase_a"_a, "angle_B"_a, "phase_b"_a,
            "Apply an RF hard pulse to an two-pools exchange model")
        .def(
            "apply_pulse", 
            static_cast<
                    void(Base::*)(Quantity const &, Quantity const &, Real)
                >(&Base::apply_pulse),
            "angle_a"_a, "phase_a"_a, "saturation"_a,
            "Apply an RF pulse to an two-pools magnetization transfer model")
        .def(
            "relaxation", &Base::relaxation, "duration"_a,
            "Simulate the relaxation during given duration")
        .def(
            "off_resonance", &Base::off_resonance, 
            "duration"_a,
            "Simulate field- and species-related off-resonance effects during "
                "given duration with given frequency offset")
        .def("__len__", &Base::size, "Number of states of the model")
    ;
}
