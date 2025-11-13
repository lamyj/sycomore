#include <sstream>

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityArray.h"
#include "sycomore/QuantityTensor.h"
#include "sycomore/QuantityTensorFixed.h"
#include "sycomore/Quantity.h"

#include "Quantity.h"

namespace
{

template<typename T>
std::vector<std::size_t> normalize_index(
    T const & magnitude, std::vector<ssize_t> const & i)
{
    if(magnitude.dimension() != i.size())
    {
        throw std::out_of_range(
            "Number of arguments (" + std::to_string(i.size())
                + ") does not match the number of dimensions ("
                + std::to_string(magnitude.dimension()) + ")");
    }
    
    // Convert signed values (counting from end) to unsigned values
    std::vector<std::size_t> signed_i(i.size());
    for(std::size_t d=0; d != i.size(); ++d)
    {
        if(i[d] < 0)
        {
            signed_i[d] = magnitude.shape()[d] + i[d];
        }
        else
        {
            signed_i[d] = i[d];
        }
    }
    
    xt::check_element_index(magnitude.shape(), signed_i.begin(), signed_i.end());
    
    return signed_i;
}

template<typename T>
sycomore::Quantity getitem(T const & q, std::vector<ssize_t> const & i)
{
    return q[normalize_index(q.magnitude, i)];
}

template<typename T>
sycomore::Quantity getitem(T const & q, ssize_t i)
{
    return getitem(q, std::vector<ssize_t>{i});
}

template<typename T>
sycomore::Quantity const & setitem(
    T & l, std::vector<ssize_t> const & i, sycomore::Quantity const & r)
{
    l[normalize_index(l.magnitude, i)] = r;
    return r;
}

template<typename T>
sycomore::Quantity const & setitem(T & l, ssize_t i, sycomore::Quantity const & r)
{
    return setitem(l, std::vector<ssize_t>{i}, r);
}

}

template<typename T>
pybind11::class_<T>
wrap_quantity(pybind11::module & m, std::string const & name)
{
    using namespace pybind11;
    using namespace sycomore;
    
    using Container = typename T::Container;
    
    auto _class = class_<T>(m, name.c_str(), "Quantity in the SI system.")
        .def(init<>())
        .def(
            init<Container, Dimensions>(),
            "magnitude"_a, "dimensions"_a=Dimensionless)
        .def_readwrite(
            "magnitude", &T::magnitude, 
            "The magnitude of the quantity, in SI units.")
        .def_readwrite(
            "dimensions", &T::dimensions, "Dimensions of the quantity")
        .def(
            "convert_to", &T::template convert_to<Quantity>, 
            "Return the scalar value of the quantity converted to the given "
            "unit.")
        .def(self == self, "Test whether magnitudes and dimensions are equal")
        .def(
            self == Container(),
            "Test whether magnitudes and dimensions are equal")
        .def(
            Container() == self,
            "Test whether magnitudes and dimensions are equal")
        .def(self != self, "Test whether magnitudes or dimensions differ")
        .def(self != Container(), "Test whether magnitudes or dimensions differ")
        .def(Container() != self, "Test whether magnitudes or dimensions differ")
        .def(self += self, "In-place addition of a compatible quantity")
        .def(self += Container(), "In-place addition of a compatible quantity")
        .def(self -= self, "In-place subtraction of a compatible quantity")
        .def(self -= Container(), "In-place subtraction of a compatible quantity")
        .def(self *= self, "In-place multiplication")
        .def(self *= Container(), "In-place multiplication")
        .def(self /= self, "In-place division")
        .def(self /= Container(), "In-place division")
        .def(
            "__ifloordiv__",
            [](T & l, T const & r) { l = std::floor(l/r); return l; },
            is_operator())
        .def(self %= Container(), "In-place floating-point modulo")
        .def(self %= self, "In-place floating-point modulo")
        .def(+self, "Identity operator")
        .def(-self, "Return a quantity with the opposite magnitude")
        .def(self + self, "Addition of compatible quantities")
        .def("__add__", sycomore::operator+<T, Container>, is_operator())
        .def(
            "__radd__",
            [](T const & r, Container const & l) { return l + r; },
            is_operator())
        .def(self - self, "Subtraction of compatible quantities")
        .def("__sub__", sycomore::operator-<T, Container>, is_operator())
        .def(
            "__rsub__",
            [](T const & r, Container const & l) { return l - r; },
            is_operator())
        .def(self * self, "Multiplication")
        .def("__mul__", sycomore::operator*<T, Container>, is_operator())
        .def(
            "__rmul__",
            [](T const & r, Container const & l) { return l * r; },
            is_operator())
        .def(self / self, "Division")
        .def("__truediv__", sycomore::operator/<T, Container>, is_operator())
        .def(
            "__rtruediv__",
            [](T const & r, Container const & l) { return l / r; },
            is_operator())
        .def(
            "__floordiv__",
            [](T const & l, T const & r) { return std::floor(l/r); },
            is_operator())
        .def(
            "__floordiv__",
            [](T const & l, Container const & r) { return std::floor(l/r); },
            is_operator())
        .def(
            "__rfloordiv__",
            [](T const & r, Container const & l) { return std::floor(l/r); },
            is_operator())
        .def("__mod__", sycomore::fmod<T, T>, is_operator())
        .def("__mod__", sycomore::fmod<T, Container>, is_operator())
        .def(
            "__divmod__",
            [](object const & l, object const & r) {
                return make_tuple(l.attr("__floordiv__")(r), l.attr("__mod__")(r));
        })
        .def(
            "__abs__", overload_cast<T const &>(std::abs<T>),
            "Return a quantity with the absolute value of the magnitude")
        .def(
            "__pow__", overload_cast<T const &, double>(std::pow<T>),
            "Raise a quantity to a power")
        .def(
            "__round__", overload_cast<T const &>(std::round<T>),
            "Round the magnitude of a quantity")
        .def(
            "__trunc__", overload_cast<T const &>(std::trunc<T>),
            "Truncate the magnitude of a quantity")
        .def(
            "__floor__", overload_cast<T const &>(std::floor<T>),
            "Quantity with the largest integer magnitude not greater than the "
                "magnitude")
        .def(
            "__ceil__", overload_cast<T const &>(std::ceil<T>),
            "Quantity with the smallest integer magnitude not less than the "
                "magnitude")
        .def(
            "__repr__",
            [](T const & d) {
                std::ostringstream s;
                s << d;
                return s.str();
            },
            "String representation of a quantity")
        .def(hash(self))
        .def(pickle(
            [](T const & q) {
                return pybind11::make_tuple(
                    q.magnitude, 
                    q.dimensions.length, q.dimensions.mass, q.dimensions.time,
                    q.dimensions.electric_current, 
                    q.dimensions.thermodynamic_temperature,
                    q.dimensions.amount_of_substance, 
                    q.dimensions.luminous_intensity);
            },
            [](pybind11::tuple t) {
                if(t.size() != 8)
                {
                    throw std::runtime_error("Invalid state!");
                }
                Dimensions dimensions(
                    t[1].cast<double>(), t[2].cast<double>(), 
                    t[3].cast<double>(), t[4].cast<double>(), 
                    t[5].cast<double>(), t[6].cast<double>(), 
                    t[7].cast<double>());
                return T(t[0].cast<typename T::Container>(), dimensions);
            }
        ));
    
    #define DIMENSIONLESS(x) (x).check_dimensions(Dimensionless)
    
    _class
        // ufuncs of numpy
        .def(
            "fmod", overload_cast<T const &, T const &>(sycomore::fmod<T, T>),
            "Floating-point modulo")
        .def(
            "fmod",
            overload_cast<T const &, Container const &>(sycomore::fmod<T, Container>),
            "Floating-point modulo")
        .def(
            "fabs", overload_cast<T const &>(std::abs<T>),
            "Return a quantity with the absolute value of the magnitude")
        .def(
            "rint", overload_cast<T const &>(std::round<T>),
            "Round the magnitude of a quantity")
        .def("exp", [](T const & x) { DIMENSIONLESS(x); return T(exp(x.magnitude));})
        .def("exp2", [](T const & x) { DIMENSIONLESS(x); return T(exp2(x.magnitude));})
        .def("log", [](T const & x) { DIMENSIONLESS(x); return T(log(x.magnitude));})
        .def("log2", [](T const & x) { DIMENSIONLESS(x); return T(log2(x.magnitude));})
        .def("log10", [](T const & x) { DIMENSIONLESS(x); return T(log10(x.magnitude));})
        .def("expm1", [](T const & x) { DIMENSIONLESS(x); return T(expm1(x.magnitude));})
        .def("log1p", [](T const & x) { DIMENSIONLESS(x); return T(log1p(x.magnitude));})
        .def("sqrt", [](T const & x) { return std::pow(x, 0.5);})
        .def("cbrt", [](T const & x) { return std::pow(x, 1./3.);})
        .def("sin", [](T const & x) { DIMENSIONLESS(x); return T(sin(x.magnitude));})
        .def("cos", [](T const & x) { DIMENSIONLESS(x); return T(cos(x.magnitude));})
        .def("tan", [](T const & x) { DIMENSIONLESS(x); return T(tan(x.magnitude));})
        .def("arcsin", [](T const & x) { DIMENSIONLESS(x); return T(asin(x.magnitude));})
        .def("arccos", [](T const & x) { DIMENSIONLESS(x); return T(acos(x.magnitude));})
        .def("arctan", [](T const & x) { DIMENSIONLESS(x); return T(atan(x.magnitude));})
        .def("arctan2", [](T const & x, T const & y) {
            DIMENSIONLESS(x); DIMENSIONLESS(y);
            return T(atan2(x.magnitude, y.magnitude));})
        .def("hypot", [](T const & x, T const & y) {
            DIMENSIONLESS(x); DIMENSIONLESS(y);
            return T(hypot(x.magnitude, y.magnitude));})
        .def("sinh", [](T const & x) { DIMENSIONLESS(x); return T(sinh(x.magnitude));})
        .def("cosh", [](T const & x) { DIMENSIONLESS(x); return T(cosh(x.magnitude));})
        .def("tanh", [](T const & x) { DIMENSIONLESS(x); return T(tanh(x.magnitude));})
        .def("arcsinh", [](T const & x) { DIMENSIONLESS(x); return T(asinh(x.magnitude));})
        .def("arccosh", [](T const & x) { DIMENSIONLESS(x); return T(acosh(x.magnitude));})
        .def("arctanh", [](T const & x) { DIMENSIONLESS(x); return T(atanh(x.magnitude));})
        .def(
            "ceil", overload_cast<T const &>(std::ceil<T>),
            "Quantity with the smallest integer magnitude not less than the "
                "magnitude")
        .def(
            "floor", overload_cast<T const &>(std::floor<T>),
            "Quantity with the largest integer magnitude not greater than the "
                "magnitude")
        .def(
            "trunc", overload_cast<T const &>(std::trunc<T>),
            "Truncate the magnitude of a quantity");
    
    #undef DIMENSIONLESS
    
    #define IN_PLACE_OPERATOR(name, op) \
        _class.def( \
            name, [](T & l, Quantity const & r) { return (l op r); }, \
            is_operator());
    #define OUT_OF_PLACE_OPERATOR(name, op) \
        _class.def( \
            name, sycomore::operator op <T, Quantity>, \
            is_operator());
    #define OUT_OF_PLACE_R_OPERATOR(name, op) \
        _class.def( \
            name, [](T const & r, Quantity const & l) { return (l op r); }, \
            is_operator());
    #define OUT_OF_PLACE_OPERATORS(name, op) \
        TQ_OPERATOR(name, op); \
        QT_OPERATOR(name, op);
    
    if(!std::is_same_v<T, Quantity>)
    {
        IN_PLACE_OPERATOR("__iadd__", +=);
        IN_PLACE_OPERATOR("__isub__", -=);
        IN_PLACE_OPERATOR("__imul__", *=);
        IN_PLACE_OPERATOR("__itruediv__", /=);
        _class.def(
            "__ifloordiv__",
            [](T & l, Quantity const & r) { l = std::floor(l/r); return l; },
            is_operator());
        IN_PLACE_OPERATOR("__imod__", %=);
        
        OUT_OF_PLACE_OPERATOR("__add__", +);
        OUT_OF_PLACE_R_OPERATOR("__radd__", +);
        OUT_OF_PLACE_OPERATOR("__sub__", -);
        OUT_OF_PLACE_R_OPERATOR("__rsub__", -);
        OUT_OF_PLACE_OPERATOR("__mul__", *);
        OUT_OF_PLACE_R_OPERATOR("__rmul__", *);
        OUT_OF_PLACE_OPERATOR("__truediv__", /);
        OUT_OF_PLACE_R_OPERATOR("__rtruediv__", /);
        _class.def(
            "__floordiv__",
            [](T const & l, Quantity const & r) { return std::floor(l/r); },
            is_operator());
        _class.def(
            "__rfloordiv__",
            [](T const & r, Quantity const & l) { return std::floor(l/r); },
            is_operator());
        _class.def(
            "__mod__",
            overload_cast<T const &, Quantity const &>(sycomore::fmod<T, Quantity>),
            is_operator());
        _class.def(
            "__mod__",
            overload_cast<T const &, double const &>(sycomore::fmod<T, double>),
            is_operator());
        _class.def(
            "fmod",
            overload_cast<T const &, Quantity const &>(sycomore::fmod<T, Quantity>),
            "Floating-point modulo");
        _class.def(
            "fmod",
            overload_cast<T const &, double const &>(sycomore::fmod<T, double>),
            "Floating-point modulo");
        
        #undef OUT_OF_PLACE_R_OPERATOR
        #undef OUT_OF_PLACE_OPERATOR
        #undef IN_PLACE_OPERATOR
    }
    
    return _class;
}

void wrap_Quantity(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    auto QuantityClass = wrap_quantity<Quantity>(m, "Quantity");
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
    
    #define WRAP_QUANTITY_CONTAINER(C) \
        auto C ## Class = wrap_quantity<C>(m, #C); \
        C ## Class \
            .def(init(&wrappers::as_quantity<C>)) \
            .def( \
                "__getitem__", \
                overload_cast<C const &, std::vector<ssize_t> const &>(getitem<C>)) \
            .def( \
                "__getitem__", \
                overload_cast<C const &, ssize_t>(getitem<C>)) \
            .def( \
                "__setitem__", \
                overload_cast<C &, std::vector<ssize_t> const &, Quantity const &>( \
                    setitem<C>)) \
            .def( \
                "__setitem__", \
                overload_cast<C &, ssize_t, Quantity const &>(setitem<C>));
            
    WRAP_QUANTITY_CONTAINER(Vector2Q);
    WRAP_QUANTITY_CONTAINER(Vector3Q);
    WRAP_QUANTITY_CONTAINER(Vector4Q);
    WRAP_QUANTITY_CONTAINER(Matrix2x2Q);
    WRAP_QUANTITY_CONTAINER(Matrix3x3Q);
    WRAP_QUANTITY_CONTAINER(ArrayQ);
}
