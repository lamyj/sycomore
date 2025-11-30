#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor/xarray.hpp>

#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityArray.h"
#include "sycomore/Quantity.h"

#include "Quantity.h"

void wrap_Quantity(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    auto QuantityClass = wrappers::wrap_quantity_class<Quantity>(m, "Quantity");
    QuantityClass
        .def(
            "__mul__",
            [](Quantity const & l, xt::xarray<double> const & r) { return l*r; },
            is_operator())
        .def(
            "__rmul__",
            [](Quantity const & r, xt::xarray<double> const & l) { return l*r; },
            is_operator())
        .def(
            "__truediv__",
            [](Quantity const & l, xt::xarray<double> const & r) { return l/r; },
            is_operator())
        .def(
            "__rtruediv__",
            [](Quantity const & r, xt::xarray<double> const & l) { return l/r; },
            is_operator())
        .def(self > self, "Compare the magnitude of two compatible quantities")
        .def(
            self > double(),
            "Compare the magnitude of two compatible quantities")
        .def(
            double() > self,
            "Compare the magnitude of two compatible quantities")
        .def(self >= self, "Compare the magnitude of two compatible quantities")
        .def(
            self >= double(),
            "Compare the magnitude of two compatible quantities")
        .def(
            double() >= self,
            "Compare the magnitude of two compatible quantities")
        .def(self < self, "Compare the magnitude of two compatible quantities")
        .def(
            self < double(),
            "Compare the magnitude of two compatible quantities")
        .def(
            double() < self,
            "Compare the magnitude of two compatible quantities")
        .def(self <= self, "Compare the magnitude of two compatible quantities")
        .def(
            self <= double(),
            "Compare the magnitude of two compatible quantities")
        .def(
            double() <= self,
            "Compare the magnitude of two compatible quantities")
        .def("__int__", [](Quantity const & q) { return int(double(q)); })
        .def(
            "__float__", [](Quantity const & q) { return double(q); },
            "Convert to a scalar");
}
