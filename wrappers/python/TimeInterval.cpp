#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/TimeInterval.h"

void wrap_TimeInterval(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<TimeInterval>(m, "TimeInterval")
        .def(
            init<Quantity, Quantity>(),
            arg("duration"), arg("gradient_moment")=0/sycomore::units::m)
        .def(init(
            [&](Quantity duration, sequence s) {
                auto sycomore_py = module::import("sycomore");
                auto Array_py = sycomore_py.attr("Array")[sycomore_py.attr("Quantity")];
                return TimeInterval(duration, Array_py(s).cast<Array<Quantity>>());
            }))
        .def_property(
            "duration",
            &TimeInterval::get_duration, &TimeInterval::set_duration)
        .def_property(
            "gradient_moment",
            &TimeInterval::get_gradient_moment,
            &TimeInterval::set_gradient_moment)
        .def(self == self)
        .def(self != self);
}
