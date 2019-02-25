#include <functional>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"

#include <iostream>

namespace std
{

template<>
struct hash<sycomore::Quantity>
{
    // Not a good hash function, but easy to test
    size_t operator()(sycomore::Quantity const & q)
    {
        auto const k1 = std::hash<double>{}(q.magnitude);
        // FIXME: k2 is in Z not in N
        auto k2 =
            (q.dimensions.length << 1)
            + (q.dimensions.mass << 2)
            + (q.dimensions.time << 3)
            + (q.dimensions.electric_current << 4)
            + (q.dimensions.thermodynamic_temperature << 5)
            + (q.dimensions.amount_of_substance << 6)
            + (q.dimensions.luminous_intensity << 7);
        if(k2 < 0)
        {
            k2 = -2*k2+1;
        }
        else
        {
            k2 = 2*k2;
        }
        // Cantor's pairing function
        return std::hash<double>{}(0.5*(k1+k2)*(k1+k2+1) + k2);
    }
};

}

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
        .def(self %= double())
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
        .def(self % double())
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
            })
        .def(hash(self));
}
