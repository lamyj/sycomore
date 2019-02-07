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
        .def(init<Real, Real>(), arg("duration")=0, arg("gradient_moment")=0)
        .def(
            init<Quantity, Quantity>(),
            arg("duration"), arg("gradient_moment")=0/sycomore::units::m)
        .def(init<Real, Array<Real>>())
        .def(init(
            [&](Real duration, sequence s) {
                auto Array_py = module::import("sycomore").attr("Array")[float_().get_type()];
                return TimeInterval(duration, Array_py(s).cast<Array<Real>>());
            }))
        .def(init(
            [&](Quantity duration, sequence s) {
                auto sycomore_py = module::import("sycomore");
                auto Array_py = sycomore_py.attr("Array")[sycomore_py.attr("Quantity")];
                return TimeInterval(duration, Array_py(s).cast<Array<Quantity>>());
            }))
        .def_readwrite("duration", &TimeInterval::duration)
        .def_readwrite("gradient_moment", &TimeInterval::gradient_moment)
        .def(self == self)
        .def(self != self);
}
