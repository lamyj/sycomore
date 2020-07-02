#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/TimeInterval.h"

void set_gradient(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<pybind11::sequence>(value))
    {
        sycomore::Array<sycomore::Quantity> array(pybind11::len(value));
        std::transform(
            value.begin(), value.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        time_interval.set_gradient(array);
    }
    else
    {
        time_interval.set_gradient(value.cast<sycomore::Quantity>());
    }
}

void set_gradient_amplitude(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<pybind11::sequence>(value))
    {
        sycomore::Array<sycomore::Quantity> array(pybind11::len(value));
        std::transform(
            value.begin(), value.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        time_interval.set_gradient_amplitude(array);
    }
    else
    {
        time_interval.set_gradient_amplitude(value.cast<sycomore::Quantity>());
    }
}

void set_gradient_area(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<pybind11::sequence>(value))
    {
        sycomore::Array<sycomore::Quantity> array(pybind11::len(value));
        std::transform(
            value.begin(), value.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        time_interval.set_gradient_area(array);
    }
    else
    {
        time_interval.set_gradient_area(value.cast<sycomore::Quantity>());
    }
}

void set_gradient_dephasing(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<pybind11::sequence>(value))
    {
        sycomore::Array<sycomore::Quantity> array(pybind11::len(value));
        std::transform(
            value.begin(), value.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        time_interval.set_gradient_dephasing(array);
    }
    else
    {
        time_interval.set_gradient_dephasing(value.cast<sycomore::Quantity>());
    }
}

void set_gradient_moment(
    sycomore::TimeInterval & time_interval, pybind11::object const & value)
{
    if(pybind11::isinstance<pybind11::sequence>(value))
    {
        sycomore::Array<sycomore::Quantity> array(pybind11::len(value));
        std::transform(
            value.begin(), value.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        time_interval.set_gradient_moment(array);
    }
    else
    {
        time_interval.set_gradient_moment(value.cast<sycomore::Quantity>());
    }
}

void wrap_TimeInterval(pybind11::module & m)
{
    using namespace pybind11;
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
        .def(init(
            [&](Quantity duration, sequence s) {
                auto sycomore_py = module::import("sycomore");
                auto Array_py = sycomore_py.attr("Array")[sycomore_py.attr("Quantity")];
                return TimeInterval(duration, Array_py(s).cast<Array<Quantity>>());
            }),
            arg("duration"), 
                arg("gradient")=Array<Quantity>{0*T_per_m, 0*T_per_m, 0*T_per_m},
            "Constructor, gradient may be specified as amplitude (in T/m), "
            "area (in T/m*s) or dephasing (in rad/m).")
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
        .def_property(
            "gradient_moment",
            &TimeInterval::get_gradient_moment, &set_gradient_moment,
            "Gradient moment, i.e. dephasing")
        .def(self == self)
        .def(self != self);
}
