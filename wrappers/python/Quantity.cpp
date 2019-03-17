#include <algorithm>
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
    size_t operator()(sycomore::Quantity const & q)
    {
        std::string binary_representation;
        std::copy(
            &q.magnitude, &q.magnitude+sizeof(q.magnitude),
            std::back_inserter(binary_representation));
        std::copy(
            &q.dimensions.length,
            &q.dimensions.length+sizeof(q.dimensions.length),
            std::back_inserter(binary_representation));
        std::copy(
            &q.dimensions.mass,
            &q.dimensions.mass+sizeof(q.dimensions.mass),
            std::back_inserter(binary_representation));
        std::copy(
            &q.dimensions.time,
            &q.dimensions.time+sizeof(q.dimensions.time),
            std::back_inserter(binary_representation));
        std::copy(
            &q.dimensions.electric_current,
            &q.dimensions.electric_current
                +sizeof(q.dimensions.electric_current),
            std::back_inserter(binary_representation));
        std::copy(
            &q.dimensions.thermodynamic_temperature,
            &q.dimensions.thermodynamic_temperature
                +sizeof(q.dimensions.thermodynamic_temperature),
            std::back_inserter(binary_representation));
        std::copy(
            &q.dimensions.amount_of_substance,
            &q.dimensions.amount_of_substance
                +sizeof(q.dimensions.amount_of_substance),
            std::back_inserter(binary_representation));
        std::copy(
            &q.dimensions.luminous_intensity,
            &q.dimensions.luminous_intensity
                +sizeof(q.dimensions.luminous_intensity),
            std::back_inserter(binary_representation));
        return std::hash<std::string>{}(binary_representation);
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
        .def(self % self)
        .def(self % double())
        .def("__pow__", static_cast<Quantity(*)(Quantity, double)>(std::pow))
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
