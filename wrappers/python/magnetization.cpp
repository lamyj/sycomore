#include <pybind11/complex.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <sycomore/magnetization.h>

void wrap_magnetization(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<ComplexMagnetization>(m, "ComplexMagnetization")
        .def(init<Complex, Real, Complex>())
        .def_readwrite("p", &ComplexMagnetization::p)
        .def_readwrite("z", &ComplexMagnetization::z)
        .def_readwrite("m", &ComplexMagnetization::m)
        .def(self == self)
        .def(self != self);
    m.attr("Magnetization") = m.attr("Array")[float_().get_type()];
    m.def("transversal", transversal);
}
