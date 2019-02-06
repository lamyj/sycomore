#include <string>

#include <pybind11/complex.h>
#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Grid.h"
#include "sycomore/magnetization.h"
#include "sycomore/sycomore.h"

template<typename T>
void wrap_Grid(
    pybind11::module & m, pybind11::handle python_type,
    std::string const & suffix)
{
    using namespace pybind11;
    using namespace sycomore;

    auto const name = "_Grid"+suffix;
    auto const grid = class_<Grid<T>>(m, name.c_str())
        .def(init<>())
        .def(init<Index, Shape>())
        .def(init<Index, Shape, T>())
        .def(
            "__getitem__",
            [](Grid<T> const & g, Index const & i) { return g[i]; })
        .def("__setitem__", [](Grid<T> & g, Index const & i, T v) { g[i]=v; })
        .def("dimension", &Grid<T>::dimension)
        .def(
            "origin", &Grid<T>::origin,
            return_value_policy::reference_internal)
        .def("shape", &Grid<T>::shape, return_value_policy::reference_internal)
        .def(
            "stride", &Grid<T>::stride,
            return_value_policy::reference_internal)
        .def(
            "__iter__",
            [](Grid<T> & g) {
                return make_iterator<
                    return_value_policy::reference_internal,
                    typename Grid<T>::iterator, typename Grid<T>::iterator,
                    T &
                >(g.begin(), g.end());
            },
            keep_alive<0, 1>());

    m.attr("Grid")[python_type] = grid;
}

void wrap_Grid(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    m.attr("Grid") = dict();

    wrap_Grid<Complex>(m, eval("complex"), "Complex");
    wrap_Grid<ComplexMagnetization>(
        m, m.attr("ComplexMagnetization"), "ComplexMagnetization");
}
