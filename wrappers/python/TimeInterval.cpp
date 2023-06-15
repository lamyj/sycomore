#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/TimeInterval.h"

#include "type_casters.h"

void set_gradient(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<sycomore::Quantity>(value))
    {
        time_interval.set_gradient(value.cast<sycomore::Quantity>());
    }
    else
    {
        time_interval.set_gradient(value.cast<sycomore::Vector3Q>());
    }
}

void set_gradient_amplitude(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<sycomore::Quantity>(value))
    {
        time_interval.set_gradient_amplitude(value.cast<sycomore::Quantity>());
    }
    else
    {
        time_interval.set_gradient_amplitude(value.cast<sycomore::Vector3Q>());
    }
}

void set_gradient_area(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<sycomore::Quantity>(value))
    {
        time_interval.set_gradient_area(value.cast<sycomore::Quantity>());
    }
    else
    {
        time_interval.set_gradient_area(value.cast<sycomore::Vector3Q>());
    }
}

void set_gradient_dephasing(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<sycomore::Quantity>(value))
    {
        time_interval.set_gradient_dephasing(value.cast<sycomore::Quantity>());
    }
    else
    {
        time_interval.set_gradient_dephasing(value.cast<sycomore::Vector3Q>());
    }
}

void wrap_TimeInterval(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    
    using namespace sycomore;
    
    using sycomore::units::T;
    using sycomore::units::Hz;
    
    auto const T_per_m = T/units::m;

    class_<TimeInterval>(
            m, "TimeInterval",  
            "Time interval, with or without magnetic field gradient.")
        .def(
            init<Quantity, Quantity>(),
            arg("duration"), arg("gradient")=0*T_per_m,
            "Constructor, gradient may be specified as amplitude (in T/m), "
            "area (in T/m*s) or dephasing (in rad/m).")
        .def(
            init<Quantity, Vector3Q>(),
            arg("duration"), arg("gradient"),
            "Constructor, gradient may be specified as amplitude (in T/m), "
            "area (in T/m*s) or dephasing (in rad/m).")
        .def_static(
            "shortest",
            static_cast<TimeInterval(*)(Quantity const &, Quantity const &)>(
                &TimeInterval::shortest),
            "k"_a, "G_max"_a,
            "Shortest possible time interval given 3D gradient area (T/m*s) "
            "or dephasing (rad/m) and maximum gradient amplitude")
        .def_static(
            "shortest",
            static_cast<TimeInterval(*)(Vector3Q const &, Quantity const &)>(
                &TimeInterval::shortest),
            "k"_a, "G_max"_a,
            "Shortest possible time interval given 3D gradient area (T/m*s) "
            "or dephasing (rad/m) and maximum gradient amplitude")
        .def_property(
            "duration",
            &TimeInterval::get_duration, &TimeInterval::set_duration, 
            "Duration")
        .def(
            "set_gradient", &set_gradient, 
            "Set gradient amplitude (in T/m), area (in T/m*s) or dephasing "
            "(in rad/m).")
        .def_property(
            "gradient_amplitude",
            &TimeInterval::get_gradient_amplitude, &set_gradient_amplitude,
            "Gradient amplitude.")
        .def_property(
            "gradient_area",
            &TimeInterval::get_gradient_area, &set_gradient_area,
            "Gradient area.")
        .def_property(
            "gradient_dephasing",
            &TimeInterval::get_gradient_dephasing, &set_gradient_dephasing,
            "Gradient dephasing")
        .def(self == self)
        .def(self != self);
}
