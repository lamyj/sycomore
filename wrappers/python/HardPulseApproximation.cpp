#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sycomore/HardPulseApproximation.h>

void wrap_HardPulseApproximation(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<HardPulseApproximation>(m, "HardPulseApproximation")
        .def(init<
            Pulse, std::vector<Quantity>, HardPulseApproximation::Envelope,
            std::string>())
        .def(init<
            Pulse, std::vector<Quantity>, HardPulseApproximation::Envelope,
            Quantity, Quantity, std::string>())
        .def("get_pulses", &HardPulseApproximation::get_pulses)
        .def("get_time_interval", &HardPulseApproximation::get_time_interval)
        .def("get_name", &HardPulseApproximation::get_name)
        .def("get_gradient_moment", &HardPulseApproximation::get_gradient_moment)
        .def("set_phase", &HardPulseApproximation::set_phase);

    m.def("sinc_envelope", sinc_envelope);
}
