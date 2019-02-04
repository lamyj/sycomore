#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Array.h"
#include "sycomore/sycomore.h"

template<typename T>
void wrap_Array(pybind11::module & m, pybind11::object python_type, std::string const & suffix)
{
    using namespace pybind11;
    using namespace sycomore;

    auto const name = "_Array"+suffix;
    auto const array = class_<Array<T>>(m, name.c_str(), buffer_protocol())
        .def(init<>())
        .def("__init__", [](Array<T> & array, buffer b) {
            buffer_info info = b.request();
            if(info.format != format_descriptor<T>::format())
            {
                throw std::runtime_error("Incompatible format");
            }
            if(info.ndim != 1)
            {
                throw std::runtime_error("Incompatible dimension");
            }
            new (&array) Array<T>(info.shape[0]);
            T * data = static_cast<T *>(info.ptr);
            std::copy(data, data+array.size(), array.begin());
        })
        .def("__init__", [](Array<T> & array, sequence s) {
            new (&array) Array<T>(s.size());
            for(size_t i=0; i<array.size(); ++i)
            {
                array[i] = s[i].cast<T>();
            }
        })
        .def_buffer([](Array<T> & array) -> buffer_info {
            return buffer_info(
                array.begin(), sizeof(T), format_descriptor<T>::format(),
                1, { array.size() }, { sizeof(T) });
        })
        .def("size", &Array<T>::size)
        .def("__len__", &Array<T>::size)
        .def("empty", &Array<T>::empty)
        .def("__bool__", [](Array<T> const & a) { return !a.empty(); })
        .def("__getitem__", [](Array<T> const & a, size_t i) { return a[i]; })
        .def("__setitem__", [](Array<T> & a, size_t i, T v) { a[i]=v; })
        .def(self == self)
        .def(self != self)
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

    m.attr("Array")[python_type.get_type()] = array;
}

void wrap_Array(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    m.attr("Array") = dict();

    wrap_Array<Real>(m, float_(), "Real");
    wrap_Array<Quantity>(m, m.attr("Quantity"), "Quantity");
}
