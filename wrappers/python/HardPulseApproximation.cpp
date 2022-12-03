#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sycomore/HardPulseApproximation.h>

void wrap_HardPulseApproximation(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    using namespace sycomore;

    class_<HardPulseApproximation>(m, "HardPulseApproximation")
        .def(
            init<
                Pulse, std::vector<Quantity>, HardPulseApproximation::Envelope,
                std::string>(),
            "model"_a, "support"_a, "envelope"_a, "name"_a="")
        .def(
            init<
                Pulse, std::vector<Quantity>, HardPulseApproximation::Envelope,
                Quantity, Quantity, std::string>(),
            "model"_a, "support"_a, "envelope"_a, "bandwidth"_a,
            "slice_thickness"_a, "name"_a="")
        .def_property_readonly(
            "pulses", &HardPulseApproximation::get_pulses)
        .def_property_readonly(
            "time_interval", &HardPulseApproximation::get_time_interval)
        .def_property_readonly(
            "name", &HardPulseApproximation::get_name)
        .def_property_readonly(
            "gradient_moment", &HardPulseApproximation::get_gradient_moment)
        .def("set_phase", &HardPulseApproximation::set_phase, "phase"_a);

    m.def("sinc_envelope", sinc_envelope);
}
