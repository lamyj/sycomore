#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <xtensor-python/pytensor.hpp>

#include "sycomore/isochromat/Model.h"

#include "../type_casters.h"

void wrap_isochromat_Model(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    using namespace sycomore;
    using namespace sycomore::isochromat;

    class_<Model>(m, "Model")
        .def(
            init<
                Quantity const &, Quantity const &, TensorR<1> const &,
                TensorQ<2> const &>(),
            "T1"_a, "T2"_a, "M0"_a, "positions"_a)
        .def(
            init<
                TensorQ<1> const &, TensorQ<1> const &, TensorR<2> const &,
                TensorQ<2> const &>(),
            "T1"_a, "T2"_a, "M0"_a, "positions"_a)
        .def(
            "build_pulse",
            overload_cast<Quantity const &, Quantity const &>(
                &Model::build_pulse, const_),
            "angle"_a, "phase"_a)
        .def(
            "build_pulse",
            overload_cast<TensorQ<1> const &, TensorQ<1> const &>(
                &Model::build_pulse, const_),
            "angle"_a, "phase"_a)
        .def(
            "build_time_interval",
            overload_cast<
                    Quantity const &, Quantity const &, TensorQ<1> const &>(
                &Model::build_time_interval, const_),
            "duration"_a, "delta_omega"_a, "gradient"_a)
        .def(
            "build_time_interval",
            overload_cast<
                    Quantity const &, TensorQ<1> const &, TensorQ<2> const &>(
                &Model::build_time_interval, const_),
            "duration"_a, "delta_omega"_a, "gradient"_a)
        .def("build_relaxation", &Model::build_relaxation, "duration"_a)
        .def(
            "build_phase_accumulation",
            overload_cast<Quantity const &>(
                &Model::build_phase_accumulation, const_),
            "angle"_a)
        .def(
            "build_phase_accumulation",
            overload_cast<TensorQ<1> const &>(
                &Model::build_phase_accumulation, const_),
            "angle"_a)
        .def("apply", &Model::apply, "operator"_a)
        .def_property_readonly("T1", &Model::T1)
        .def_property_readonly("T2", &Model::T2)
        .def_property_readonly("M0", &Model::M0)
        .def_property_readonly("magnetization", &Model::magnetization)
        .def_property_readonly("positions", &Model::positions);
}
