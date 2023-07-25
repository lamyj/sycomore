#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <xtensor-python/pytensor.hpp>

#include "sycomore/isochromat/Operator.h"

#include "../type_casters.h"

void wrap_isochromat_Operator(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    using namespace sycomore;
    using namespace sycomore::isochromat;

    class_<Operator>(
            m, "Operator", 
            "Isochromat simulation operator, i.e. an array of 4×4 matrices")
        .def(init<>())
        .def(init<TensorR<3> const &>(), "data"_a)
        .def(
            self *= self,
            "In-place chaining of operators: with right applied first.\n"
            "If self and/or right has shape 1×4×4, it is broadcast to match "
            "the other operand. Otherwise, both operand must have shape n×4×4.")
        .def(
            "preMultiply", &Operator::preMultiply, "left"_a,
            "In-place chaining of operators, with self applied first.\n"
            "If self and/or left has shape 1×4×4, it is broadcast to match "
            "the other operand. Otherwise, both operand must have shape n×4×4.")
        .def_property_readonly(
            "array", [](Operator const & o){
                return xt::pytensor<Real,3>(o.array());
            },
            "Numeric representation of the operator")
        .def(
            self * self,
            "Operator chaining, representing right followed by left");
}
