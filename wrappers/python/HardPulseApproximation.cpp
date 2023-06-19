#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sycomore/HardPulseApproximation.h>

#include "type_casters.h"

void wrap_HardPulseApproximation(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    using namespace sycomore;

    class_<HardPulseApproximation>(m, "HardPulseApproximation")
        .def(
            init<
                Pulse, std::vector<Quantity>,
                HardPulseApproximation::Envelope>(),
            "model"_a, "support"_a, "envelope"_a)
        .def_property_readonly(
            "pulses", &HardPulseApproximation::pulses)
        .def_property_readonly(
            "duration", &HardPulseApproximation::duration)
        .def_property(
            "phase",
            &HardPulseApproximation::set_phase,
            &HardPulseApproximation::set_phase);

    m.def("apodized_sinc_envelope", apodized_sinc_envelope);
    m.def("hamming_sinc_envelope", hamming_sinc_envelope);
    m.def("hanning_sinc_envelope", hanning_sinc_envelope);
    m.def("sinc_envelope", sinc_envelope);
}
