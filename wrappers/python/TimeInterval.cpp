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
        .def(init<Real, Real>())
        .def(init<Quantity, Quantity>())
        .def(
            "__init__",
            [&](TimeInterval & self, Real duration, sequence s) {
                auto Array_py = module::import("sycomore").attr("Array")[float_().get_type()];
                new (&self) TimeInterval(duration, Array_py(s).cast<Array<Real>>());
            })
        .def(
            "__init__",
            [&](TimeInterval & self, Quantity duration, sequence s) {
                auto sycomore_py = module::import("sycomore");
                auto Array_py = sycomore_py.attr("Array")[sycomore_py.attr("Quantity").get_type()];
                new (&self) TimeInterval(duration, Array_py(s).cast<Array<Quantity>>());
            })
        .def_readwrite("duration", &TimeInterval::duration)
        .def_readwrite("gradient_moment", &TimeInterval::gradient_moment)
        .def(self == self)
        .def(self != self);
}
