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

    class_<TimeInterval>(m, "TimeInterval")
        .def(
            init<Quantity, Quantity>(),
            arg("duration"), arg("gradient")=0*units::T/units::m)
        .def(init(
            [&](Quantity duration, sequence s) {
                auto sycomore_py = module::import("sycomore");
                auto Array_py = sycomore_py.attr("Array")[sycomore_py.attr("Quantity")];
                return TimeInterval(duration, Array_py(s).cast<Array<Quantity>>());
            }))
        .def_property(
            "duration",
            &TimeInterval::get_duration, &TimeInterval::set_duration)
        .def("set_gradient", &set_gradient)
        .def_property(
            "gradient_amplitude",
            &TimeInterval::get_gradient_amplitude, &set_gradient_amplitude)
        .def_property(
            "gradient_area",
            &TimeInterval::get_gradient_area, &set_gradient_area)
        .def_property(
            "gradient_dephasing",
            &TimeInterval::get_gradient_dephasing, &set_gradient_dephasing)
        .def_property(
            "gradient_moment",
            &TimeInterval::get_gradient_moment, &set_gradient_moment)
        .def(self == self)
        .def(self != self);
}
