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

using UnaryFunction = double (*)(double);
template<UnaryFunction function>
sycomore::Quantity scalar_unary_function(sycomore::Quantity const & q)
{
    if(q.dimensions != sycomore::Dimensions())
    {
        throw std::runtime_error("Defined only for scalars");
    }
    return sycomore::Quantity(function(q.magnitude), q.dimensions);
}

using BinaryFunction = double (*)(double, double);
template<BinaryFunction function>
sycomore::Quantity scalar_binary_function(
    sycomore::Quantity const & q1, sycomore::Quantity const & q2)
{
    if(
        q1.dimensions != sycomore::Dimensions() 
        || q2.dimensions != sycomore::Dimensions())
    {
        throw std::runtime_error("Defined only for scalars");
    }
    return sycomore::Quantity(
        function(q1.magnitude, q2.magnitude), sycomore::Dimensions());
}

}

void wrap_Quantity(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    auto QuantityClass = class_<Quantity>(m, "Quantity")
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
        .def(self == double())
        .def(double() == self)
        .def(self != self)
        .def(self != double())
        .def(double() != self)
        .def(self += self)
        .def(self += double())
        .def(self -= self)
        .def(self -= double())
        .def(self *= self)
        .def(self *= double())
        .def(self /= self)
        .def(self /= double())
        .def(self %= double())
        .def(self %= self)
        .def(+self)
        .def(-self)
        .def(self + self)
        .def(self + double())
        .def(double() + self)
        .def(self - self)
        .def(self - double())
        .def(double() - self)
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
        .def(self > double())
        .def(double() > self)
        .def(self >= self)
        .def(self >= double())
        .def(double() >= self)
        .def(self < self)
        .def(self < double())
        .def(double() < self)
        .def(self <= self)
        .def(self <= double())
        .def(double() <= self)
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
        ))
        .def("__int__", [](Quantity const & q) { return int(double(q)); })
        .def("__float__", [](Quantity const & q) { return double(q); })
        // ufuncs of numpy
        .def("fmod", [](Quantity const & l, Quantity const & r) { return l%r; })
        .def("fmod", [](Quantity const & l, double r) { return l%r; })
        .def("fabs", static_cast<Quantity(*)(Quantity)>(std::abs))
        .def("rint", static_cast<Quantity(*)(Quantity)>(std::round))
        .def("exp", scalar_unary_function<std::exp>)
        .def("exp2", scalar_unary_function<std::exp2>)
        .def("log", scalar_unary_function<std::log>)
        .def("log2", scalar_unary_function<std::log2>)
        .def("log10", scalar_unary_function<std::log10>)
        .def("expm1", scalar_unary_function<std::expm1>)
        .def("log1p", scalar_unary_function<std::log1p>)
        .def("sqrt", [](Quantity const & q) { return std::pow(q, 0.5); })
        .def("cbrt", [](Quantity const & q) { return std::pow(q, 1./3.); })
        .def("sin", scalar_unary_function<std::sin>)
        .def("cos", scalar_unary_function<std::cos>)
        .def("tan", scalar_unary_function<std::tan>)
        .def("arcsin", scalar_unary_function<std::asin>)
        .def("arccos", scalar_unary_function<std::acos>)
        .def("arctan", scalar_unary_function<std::atan>)
        .def("arctan2", scalar_binary_function<std::atan2>)
        .def("hypot", scalar_binary_function<std::hypot>)
        .def("sinh", scalar_unary_function<std::sinh>)
        .def("cosh", scalar_unary_function<std::cosh>)
        .def("tanh", scalar_unary_function<std::tanh>)
        .def("arcsinh", scalar_unary_function<std::asinh>)
        .def("arccosh", scalar_unary_function<std::acosh>)
        .def("arctanh", scalar_unary_function<std::atanh>)
        .def("ceil", static_cast<Quantity(*)(Quantity)>(std::ceil))
        .def("floor", static_cast<Quantity(*)(Quantity)>(std::floor))
        .def("trunc", static_cast<Quantity(*)(Quantity)>(std::trunc))
    ;
}
