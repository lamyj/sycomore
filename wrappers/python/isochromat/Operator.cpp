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
        .def(self *= self)
        .def("preMultiply", &Operator::preMultiply)
        .def_property_readonly(
            "array", [](Operator const & o){
                return xt::pytensor<Real,3>(o.array());
            })
        .def(self * self);
}
