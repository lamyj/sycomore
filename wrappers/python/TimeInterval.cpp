#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/TimeInterval.h"

#include "Quantity.h"

sycomore::TimeInterval constructor(
    sycomore::Quantity const & duration, pybind11::object gradient)
{
    static auto const quantity_type =
        pybind11::module::import("sycomore").attr("Quantity");
    if(pybind11::type::of(gradient).is(quantity_type))
    {
        return {duration, gradient.cast<sycomore::Quantity>()};
    }
    else
    {
        return {
            duration,
            sycomore::wrappers::as_quantity<sycomore::Vector3Q>(gradient)};
    }
}

sycomore::TimeInterval shortest(
    pybind11::object k, sycomore::Quantity const & G_max)
{
    static auto const quantity_type =
        pybind11::module::import("sycomore").attr("Quantity");
    if(pybind11::type::of(k).is(quantity_type))
    {
        return sycomore::TimeInterval::shortest(
            k.cast<sycomore::Quantity>(), G_max);
    }
    else
    {
        return sycomore::TimeInterval::shortest(
            sycomore::wrappers::as_quantity<sycomore::Vector3Q>(k), G_max);
    }
}

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
            init(&constructor),
            arg("duration"), arg("gradient")=0*T_per_m,
            "Constructor, gradient may be specified as amplitude (in T/m), "
            "area (in T/m*s) or dephasing (in rad/m).")
        .def_static(
            "shortest", &shortest,
            "k"_a, "G_max"_a,
            "Shortest possible time interval given 3D gradient area (T/m*s) "
            "or dephasing (rad/m) and maximum gradient amplitude")
        .def_property(
            "duration", &TimeInterval::duration, &TimeInterval::set_duration, 
            "Duration")
        .def(
            "set_gradient", &set_gradient, 
            "Set gradient amplitude (in T/m), area (in T/m*s) or dephasing "
            "(in rad/m).")
        .def_property(
            "gradient_amplitude",
            &TimeInterval::gradient_amplitude, &set_gradient,
            "Gradient amplitude.")
        .def_property(
            "gradient_area", &TimeInterval::gradient_area, &set_gradient,
            "Gradient area.")
        .def_property(
            "gradient_dephasing",
            &TimeInterval::gradient_dephasing, &set_gradient,
            "Gradient dephasing")
        .def(self == self)
        .def(self != self);
}
