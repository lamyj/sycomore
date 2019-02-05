#include <pybind11/pybind11.h>

#include "sycomore/GridScanner.h"

void wrap_GridScanner(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<GridScanner>(m, "GridScanner")
        .def(init<Index, Shape>())
        .def(init<Index, Shape, Index, Shape>())
        .def(
            "__iter__",
            [](GridScanner & g) {
                return make_iterator<
                    return_value_policy::reference_internal,
                    GridScanner, GridScanner, GridScanner::value_type const &
                >(g.begin(), g.end());
            },
            keep_alive<0, 1>());
}
