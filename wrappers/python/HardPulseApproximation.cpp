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

    class_<HardPulseApproximation>(
            m, "HardPulseApproximation",
            "Small tip angle approximation of a shaped pulse")
        .def(
            init<
                Pulse, std::vector<Quantity>,
                HardPulseApproximation::Envelope>(),
            "model"_a, "support"_a, "envelope"_a,
            "Create a shaped pulse with a flip angle equivalent to given "
                "hard pulse")
        .def_property_readonly(
            "pulses", &HardPulseApproximation::pulses,
            "Return the hard pulses approximating the shaped pulse")
        .def_property_readonly(
            "duration", &HardPulseApproximation::duration,
            "Return the duration of a hard pulse")
        .def_property(
            "phase",
            &HardPulseApproximation::set_phase,
            &HardPulseApproximation::set_phase,
            "Return the phase of the shaped pulse");

    m.def(
        "apodized_sinc_envelope", apodized_sinc_envelope,
        "Create an apodized sinc envelope");
    m.def(
        "hann_sinc_envelope", hann_sinc_envelope,
        "Create an Hann-apodized sinc envelope");
    m.def(
        "hamming_sinc_envelope", hamming_sinc_envelope,
        "Create an Hamming-apodized sinc envelope");
    m.def("sinc_envelope", sinc_envelope, "Create a sinc envelope");
}
