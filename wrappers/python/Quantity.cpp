#include <algorithm>
#include <functional>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"

#include <iostream>

namespace
{

sycomore::Quantity floordiv(
    sycomore::Quantity l, sycomore::Quantity const & r)
{
    l /= r;
    l.magnitude = std::floor(l.magnitude);
    return l;
}

sycomore::Quantity floordiv(sycomore::Quantity q, double s)
{
    q /= s;
    q.magnitude = std::floor(q.magnitude);
    return q;
}

sycomore::Quantity rfloordiv(sycomore::Quantity const & q, double s)
{
    auto r = s/q;
    r.magnitude = std::floor(r.magnitude);
    return r;
}

pybind11::object divmod(pybind11::object const & l, pybind11::object const & r)
{
    return pybind11::make_tuple(l.attr("__floordiv__")(r), l.attr("__mod__")(r));
}

}

void wrap_Quantity(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<Quantity>(m, "Quantity")
        .def(init<double, Dimensions>())
        .def_readwrite(
            "magnitude", &Quantity::magnitude, 
            "The magnitude of the quantity, in SI units.")
        .def_readwrite("dimensions", &Quantity::dimensions)
        .def(
            "convert_to", &Quantity::convert_to, 
            "Return the scalar value of the quantity converted to the given "
            "unit.")
        .def(self == self)
        .def(self != self)
        .def(self += self)
        .def(self -= self)
        .def(self *= self)
        .def(self *= double())
        .def(self /= self)
        .def(self /= double())
        .def(self %= double())
        .def(self %= self)
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
        .def("__floordiv__", static_cast<Quantity(*)(Quantity, Quantity const &)>(floordiv))
        .def("__floordiv__", static_cast<Quantity(*)(Quantity, double)>(floordiv))
        .def("__rfloordiv__", static_cast<Quantity(*)(Quantity const &, double)>(rfloordiv))
        .def(self % self)
        .def(self % double())
        .def("__divmod__", divmod)
        .def("__abs__", static_cast<Quantity(*)(Quantity)>(std::abs))
        .def("__pow__", static_cast<Quantity(*)(Quantity, double)>(std::pow))
        .def("__round__", static_cast<Quantity(*)(Quantity)>(std::round))
        .def("__trunc__", static_cast<Quantity(*)(Quantity)>(std::trunc))
        .def("__floor__", static_cast<Quantity(*)(Quantity)>(std::floor))
        .def("__ceil__", static_cast<Quantity(*)(Quantity)>(std::ceil))
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
        .def(hash(self))
        .def(pickle(
            [](Quantity const & q) {
                return make_tuple(
                    q.magnitude, 
                    q.dimensions.length, q.dimensions.mass, q.dimensions.time,
                    q.dimensions.electric_current, 
                    q.dimensions.thermodynamic_temperature,
                    q.dimensions.amount_of_substance, 
                    q.dimensions.luminous_intensity);
            },
            [](tuple t) {
                if(t.size() != 8)
                {
                    throw std::runtime_error("Invalid state!");
                }
                Dimensions dimensions(
                    t[1].cast<double>(), t[2].cast<double>(), 
                    t[3].cast<double>(), t[4].cast<double>(), 
                    t[5].cast<double>(), t[6].cast<double>(), 
                    t[7].cast<double>());
                return Quantity(t[0].cast<double>(), dimensions);
            }
        ));
}
