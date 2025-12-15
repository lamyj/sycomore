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

namespace sycomore
{
namespace wrappers
{

std::size_t normalize_index(std::size_t shape, ssize_t i)
{
    std::size_t unsigned_i;
    if(i < 0)
    {
        unsigned_i = shape + i;
    }
    else
    {
        unsigned_i = i;
    }
    
    if(unsigned_i >= shape)
    {
        throw std::out_of_range(
            "index " + std::to_string(i) + " is out of bounds with size "
            + std::to_string(shape));
    }
    
    return unsigned_i;
}

}
}

void wrap_Quantity(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    auto QuantityClass = wrappers::wrap_quantity_class<Quantity>(m, "Quantity");
    QuantityClass
        .def(self * xt::xarray<double>()).def(xt::xarray<double>() * self)
        .def(self / xt::xarray<double>()).def(xt::xarray<double>() / self)
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
