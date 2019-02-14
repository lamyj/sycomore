#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"

void wrap_Quantity(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<Quantity>(m, "Quantity")
        .def(init<double, Dimensions>())
        .def_readwrite("magnitude", &Quantity::magnitude)
        .def_readwrite("dimensions", &Quantity::dimensions)
        .def("convert_to", &Quantity::convert_to)
        .def(self == self)
        .def(self != self)
        .def(self += self)
        .def(self -= self)
        .def(self *= self)
        .def(self *= double())
        .def(self /= self)
        .def(self /= double())
        .def(+self)
        .def(-self)
        .def(self + self)
        .def(self - self)
        .def(self * self)
        .def(self * double())
        .def(double() * self)
        .def(self / self)
        .def(self / double())
        .def(double() / self)
        .def(self > self)
        .def(self >= self)
        .def(self < self)
        .def(self <= self)
        .def(
            "__repr__",
            [](Quantity const & d) {
                std::ostringstream s;
                s << d;
                return s.str();
            });
}
