#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Array.h"
#include "sycomore/sycomore.h"

template<typename T>
pybind11::class_<sycomore::Array<T>>
wrap_Array(
    pybind11::module & m, pybind11::handle python_type,
    std::string const & suffix)
{
    using namespace pybind11;
    using namespace sycomore;

    auto const name = "_Array"+suffix;
    auto array = class_<Array<T>>(m, name.c_str(), buffer_protocol())
        .def(init<>())
        .def(init(
            [](sequence s) {
                Array<T> array(s.size());
                for(std::size_t i=0; i<array.size(); ++i)
                {
                    array[i] = s[i].cast<T>();
                }
                return array;
            }))
        .def(init(
            [](args s) {
                Array<T> array(s.size());
                for(std::size_t i=0; i<array.size(); ++i)
                {
                    array[i] = s[i].cast<T>();
                }
                return array;
            }))
        .def_buffer([](Array<T> & array) -> buffer_info {
            return buffer_info(
                array.begin(), sizeof(T), format_descriptor<T>::format(),
                1, { array.size() }, { sizeof(T) });
        })
        .def("size", &Array<T>::size)
        .def("__len__", &Array<T>::size)
        .def("empty", &Array<T>::empty)
        .def("__bool__", [](Array<T> const & a) { return !a.empty(); })
        .def("__getitem__", [](Array<T> const & a, std::size_t i) { return a[i]; })
        .def("__setitem__", [](Array<T> & a, std::size_t i, T v) { a[i]=v; })
        .def(self == self)
        .def(self != self)
        .def(-self)
        .def(self+self)
        .def(self-self)
        .def(self*T())
        .def(self/T())
        .def(
            "__iter__",
            [](Array<T> & a) {
                return make_iterator<
                    return_value_policy::reference_internal,
                    typename Array<T>::iterator, typename Array<T>::iterator,
                    T &
                >(a.begin(), a.end());
            },
            keep_alive<0, 1>())
        .def(
            "__repr__",
            [](Array<T> const & a) {
                std::ostringstream stream;
                stream << a;
                return stream.str();
            });

    m.attr("Array")[python_type] = array;

    return array;
}

void wrap_Array(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    m.attr("Array") = dict();

    auto int_ = module::import("numpy").attr("int32");
    auto uint_ = module::import("numpy").attr("uint32");

    wrap_Array<Real>(m, float_().get_type(), "Real");
    wrap_Array<int32_t>(m, int_, "Int");
    wrap_Array<uint32_t>(m, uint_, "UInt");
    wrap_Array<Quantity>(m, m.attr("Quantity"), "Quantity")
        .def(self *= double())
        .def(self /= double())
        .def(self * double())
        .def(double() * self)
        .def(self / double());

    m.attr("Index") = m.attr("Array")[int_];
    m.attr("Shape") = m.attr("Array")[uint_];
    m.attr("Stride") = m.attr("Array")[uint_];
    m.attr("Point") = m.attr("Array")[m.attr("Quantity")];
}
