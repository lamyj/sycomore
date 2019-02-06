#include <pybind11/complex.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <sycomore/magnetization.h>

void wrap_magnetization(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<Magnetization>(m, "Magnetization")
        .def(init<Real, Real, Real>())
        .def_readwrite("x", &Magnetization::x)
        .def_readwrite("y", &Magnetization::y)
        .def_readwrite("z", &Magnetization::z)
        .def("transversal", &Magnetization::transversal)
        .def(self == self)
        .def(self != self);

    class_<ComplexMagnetization>(m, "ComplexMagnetization")
        .def(init<Complex, Real, Complex>())
        .def_readwrite("p", &ComplexMagnetization::p)
        .def_readwrite("z", &ComplexMagnetization::z)
        .def_readwrite("m", &ComplexMagnetization::m)
        .def(self == self)
        .def(self != self);
}
