#ifndef _e0c9f20a_f96b_4c09_8459_f2be8e7df213
#define _e0c9f20a_f96b_4c09_8459_f2be8e7df213

#include <sstream>

#include "Quantity.h"

#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor/xexception.hpp>

#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"

namespace sycomore
{

namespace wrappers
{

template<typename T>
pybind11::class_<T>
wrap_quantity_class(pybind11::module & m, std::string const & name)
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
    
    /**************************************************************************/
    /************************** COMPARISON OPERATORS **************************/
    /**************************************************************************/
    _class
        .def(self == self, "Test whether magnitudes and dimensions are equal")
        .def(
            self == Container(),
            "Test whether magnitudes and dimensions are equal")
        .def(
            Container() == self,
            "Test whether magnitudes and dimensions are equal")
        .def(self != self, "Test whether magnitudes or dimensions differ")
        .def(self != Container(), "Test whether magnitudes or dimensions differ")
        .def(Container() != self, "Test whether magnitudes or dimensions differ");
    
    /**************************************************************************/
    /************************** ARITHMETIC OPERATORS **************************/
    /******************************* (IN PLACE) *******************************/
    /**************************************************************************/
    _class
        .def(self += self, "In-place addition of a compatible quantity")
        .def(self += Container(), "In-place addition of a compatible quantity")
        .def(self -= self, "In-place subtraction of a compatible quantity")
        .def(self -= Container(), "In-place subtraction of a compatible quantity")
        .def(self *= self, "In-place multiplication")
        .def(self *= Container(), "In-place multiplication")
        .def(self /= self, "In-place division")
        .def(self /= Container(), "In-place division")
        .def(self %= Container(), "In-place floating-point modulo")
        .def(self %= self, "In-place floating-point modulo")
        .def(
                "__ifloordiv__",
                [](T & l, T const & r) { l = std::floor(l/r); return l; },
                is_operator());
    
    /**************************************************************************/
    /************************** ARITHMETIC OPERATORS **************************/
    /********************** (OUT OF PLACE, HOMOGENENOUS) **********************/
    /**************************************************************************/
    _class
        .def(+self, "Identity operator")
        .def(-self, "Return a quantity with the opposite magnitude")
        .def(self + self, "Addition of compatible quantities")
        .def(self - self, "Subtraction of compatible quantities")
        .def(self * self, "Multiplication")
        .def(self / self, "Division")
        .def(
            "__floordiv__",
            [](T const & l, T const & r) { return std::floor(l/r); },
            is_operator())
        .def("__mod__", sycomore::fmod<T, T>, is_operator());
    
    /**************************************************************************/
    /************************** ARITHMETIC OPERATORS **************************/
    /********************* (OUT OF PLACE,  CLASS ON LEFT) *********************/
    /**************************************************************************/
    _class
        .def("__add__", sycomore::operator+<T, Container>, is_operator())
        .def("__sub__", sycomore::operator-<T, Container>, is_operator())
        .def("__mul__", sycomore::operator*<T, Container>, is_operator())
        .def("__truediv__", sycomore::operator/<T, Container>, is_operator())
        .def(
            "__floordiv__",
            [](T const & l, Container const & r) { return std::floor(l/r); },
            is_operator())
        .def("__mod__", sycomore::fmod<T, Container>, is_operator());
    
    /**************************************************************************/
    /************************** ARITHMETIC OPERATORS **************************/
    /********************* (OUT OF PLACE, CLASS ON RIGHT) *********************/
    /**************************************************************************/
    _class
        .def(
            "__radd__",
            [](T const & r, Container const & l) { return l + r; },
            is_operator())
        .def(
            "__rsub__",
            [](T const & r, Container const & l) { return l - r; },
            is_operator())
        .def(
            "__rmul__",
            [](T const & r, Container const & l) { return l * r; },
            is_operator())
        .def(
            "__rtruediv__",
            [](T const & r, Container const & l) { return l / r; },
            is_operator())
        .def(
            "__rfloordiv__",
            [](T const & r, Container const & l) { return std::floor(l/r); },
            is_operator())
        .def(
            "__divmod__",
            [](object const & l, object const & r) {
                return make_tuple(l.attr("__floordiv__")(r), l.attr("__mod__")(r));
        });
    
    /**************************************************************************/
    /**************** OTHER METHODS TO EMULATE NUMERIC OBJECTS ****************/
    /**************************************************************************/
    _class
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
        .def_static("__array_ufunc__", wrap_ufuncs<T>);
    
    return _class;
}

#define DIMENSIONLESS(x) (x).check_dimensions(Dimensionless, "Ufunc requires dimensionless");
#define UFUNC(expr) return cast(expr);
#define OPERATOR_UFUNC(op) \
    try \
    { \
        UFUNC(x op args[1].cast<T>()) \
    } \
    catch(std::runtime_error & e) \
    { \
        try \
        { \
            UFUNC(x op args[1].cast<Container>()) \
        } \
        catch(std::runtime_error & e) \
        { \
            UFUNC(x op args[1].cast<Quantity>()) \
        } \
    }

template<typename T>
pybind11::object
wrap_ufuncs(
    T const &, pybind11::object ufunc, std::string const & method,
    pybind11::args args, pybind11::kwargs kwargs)
{
    using namespace pybind11;
    using namespace sycomore;
    
    using Container = typename T::Container;
    
    object NotImplemented = module_::import("builtins").attr("NotImplemented");
    
    if(method != "__call__")
    {
        return NotImplemented;
    }
    
    auto const n = ufunc.attr("__name__").cast<std::string>();
    auto const x = args[0].cast<T>();
    
    if(n == "add")               {                  OPERATOR_UFUNC(+) }
    else if(n == "subtract")     {                  OPERATOR_UFUNC(-) }
    else if(n == "multiply")     {                  OPERATOR_UFUNC(*) }
    else if(n == "divide")       {                  OPERATOR_UFUNC(/) }
    else if(n == "floor_divide") {
        try
        {
            return cast(std::floor(x/args[1].cast<T>()));
        }
        catch(std::runtime_error & e)
        {
            try
            {
                return cast(std::floor(x/args[1].cast<Container>()));
            }
            catch(std::runtime_error & e)
            {
                return cast(std::floor(x/args[1].cast<Quantity>()));
            }
        }
    }
    else if(n == "absolute")     {                  UFUNC(std::abs(x)) }
    else if(n == "fabs")         {                  UFUNC(std::abs(x)) }
    else if(n == "rint")         {                  UFUNC(std::round(x)) }
    else if(n == "negative")     {                  UFUNC(-x) }
    else if(n == "positive")     {                  UFUNC(+x) }
    else if(n == "power")        {                  UFUNC(std::pow(x, args[1].cast<double>())) }
    else if(n == "remainder" || n == "fmod") {
        try
        {
            return cast(sycomore::fmod(x, args[1].cast<Quantity>()));
        }
        catch(std::runtime_error & e)
        {
            return cast(sycomore::fmod(x, args[1].cast<double>()));
        }
    }
    else if(n == "absolute")     {                  UFUNC(std::abs(x)) }
    else if(n == "exp")          { DIMENSIONLESS(x) UFUNC(T(exp(x.magnitude))) }
    else if(n == "exp2")         { DIMENSIONLESS(x) UFUNC(T(exp2(x.magnitude))) }
    else if(n == "log")          { DIMENSIONLESS(x) UFUNC(T(log(x.magnitude))) }
    else if(n == "log2")         { DIMENSIONLESS(x) UFUNC(T(log2(x.magnitude))) }
    else if(n == "log10")        { DIMENSIONLESS(x) UFUNC(T(log10(x.magnitude))) }
    else if(n == "expm1")        { DIMENSIONLESS(x) UFUNC(T(expm1(x.magnitude))) }
    else if(n == "log1p")        { DIMENSIONLESS(x) UFUNC(T(log1p(x.magnitude))) }
    else if(n == "sqrt")         {                  UFUNC(std::pow(x, 0.5)) }
    else if(n == "square")       {                  UFUNC(std::pow(x, 2.)) }
    else if(n == "cbrt")         {                  UFUNC(std::pow(x, 1./3.)) }
    else if(n == "reciprocal")   {                  UFUNC(1./x) }
    else if(n == "sin")          { DIMENSIONLESS(x) UFUNC(T(sin(x.magnitude))) }
    else if(n == "cos")          { DIMENSIONLESS(x) UFUNC(T(cos(x.magnitude))) }
    else if(n == "tan")          { DIMENSIONLESS(x) UFUNC(T(tan(x.magnitude))) }
    else if(n == "arcsin")       { DIMENSIONLESS(x) UFUNC(T(asin(x.magnitude))) }
    else if(n == "arccos")       { DIMENSIONLESS(x) UFUNC(T(acos(x.magnitude))) }
    else if(n == "arctan")       { DIMENSIONLESS(x) UFUNC(T(atan(x.magnitude))) }
    else if(n == "arctan2")      {
        DIMENSIONLESS(x)
        auto const & y = args[1].cast<T const &>(); DIMENSIONLESS(y)
        UFUNC(T(atan2(x.magnitude, y.magnitude))) }
    else if(n == "hypot")    {
        DIMENSIONLESS(x)
        auto const & y = args[1].cast<T const &>(); DIMENSIONLESS(y)
        UFUNC(T(hypot(x.magnitude, y.magnitude))) }
    else if(n == "sinh")       { DIMENSIONLESS(x) UFUNC(T(sinh(x.magnitude))) }
    else if(n == "cosh")       { DIMENSIONLESS(x) UFUNC(T(cosh(x.magnitude))) }
    else if(n == "tanh")       { DIMENSIONLESS(x) UFUNC(T(tanh(x.magnitude))) }
    else if(n == "arcsinh")    { DIMENSIONLESS(x) UFUNC(T(asinh(x.magnitude))) }
    else if(n == "arccosh")    { DIMENSIONLESS(x) UFUNC(T(acosh(x.magnitude))) }
    else if(n == "arctanh")    { DIMENSIONLESS(x) UFUNC(T(atanh(x.magnitude))) }
    else if(n == "ceil")       {                  UFUNC(std::ceil(x)) }
    else if(n == "floor")      {                  UFUNC(std::floor(x)) }
    else if(n == "trunc")      {                  UFUNC(std::trunc(x)) }
    else
    {
        print(n);
        print(args);
        print(kwargs);
        return NotImplemented;
    }
}

#undef DIMENSIONLESS
#undef UFUNC
#undef OPERATOR_UFUNC

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

template<typename T>
pybind11::class_<T>
wrap_quantity_array(pybind11::class_<T> & _class)
{
    using namespace pybind11;
    using namespace sycomore;
    
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
    
    _class
        .def(init(&wrappers::as_quantity<T>))
        .def(
            "__floordiv__",
            [](T const & l, Quantity const & r) { return std::floor(l/r); },
            is_operator())
        .def(
            "__rfloordiv__",
            [](T const & r, Quantity const & l) { return std::floor(l/r); },
            is_operator())
        .def(
            "__mod__",
            overload_cast<T const &, Quantity const &>(sycomore::fmod<T, Quantity>),
            is_operator())
        .def(
            "__mod__",
            overload_cast<T const &, double const &>(sycomore::fmod<T, double>),
            is_operator())
        .def(
            "fmod",
            overload_cast<T const &, Quantity const &>(sycomore::fmod<T, Quantity>),
            "Floating-point modulo")
        .def(
            "fmod",
            overload_cast<T const &, double const &>(sycomore::fmod<T, double>),
            "Floating-point modulo")
        .def(
            "__getitem__",
            overload_cast<T const &, std::vector<ssize_t> const &>(getitem<T>))
        .def(
            "__getitem__",
            overload_cast<T const &, ssize_t>(getitem<T>))
        .def(
            "__setitem__",
            overload_cast<T &, std::vector<ssize_t> const &, Quantity const &>(
                setitem<T>))
        .def(
            "__setitem__",
            overload_cast<T &, ssize_t, Quantity const &>(setitem<T>))
        .def(
            "__len__", [](T const & c) {
                if(c.magnitude.dimension() > 0)
                {
                    return c.shape().front();
                }
                else
                {
                    throw pybind11::type_error("len() of unsized object");
                }
            })
        .def_property_readonly(
            "shape", [](T const & c) {
                auto const shape = c.shape();
                tuple result(shape.size());
                for(std::size_t i=0; i!=shape.size(); ++i)
                {
                    PyTuple_SET_ITEM(result.ptr(), i, pybind11::cast(shape[i]).release().ptr());
                }
                return result;
            })
        .def(
            "__iter__", [](T & v) {
                return make_iterator(v.begin(), v.end());
            },
            keep_alive<0, 1>());
    
    return _class;
}

template<typename T>
T as_quantity(pybind11::array_t<pybind11::object> array)
{
    std::vector<size_t> const shape{array.shape(), array.shape()+array.ndim()};
    
    T destination(T::Container::from_shape(shape));
    auto dest_it = destination.magnitude.begin();
    
    array.resize({array.size()});
    
    if(array.size() != 0)
    {
        destination.dimensions = array.data()->cast<sycomore::Quantity>().dimensions;
    }
    
    for(auto && source: array)
    {
        auto const & q = source.cast<sycomore::Quantity>();
        destination.check_dimensions(q, "Constructor requires same dimensions");
        *dest_it = q.magnitude;
        ++dest_it;
    }
    
    array.resize(shape);
    
    return destination;
}

template<typename T>
std::vector<std::size_t>
normalize_index(T const & magnitude, std::vector<ssize_t> const & i)
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
sycomore::Quantity
getitem(T const & q, std::vector<ssize_t> const & i)
{
    return q[normalize_index(q.magnitude, i)];
}

template<typename T>
sycomore::Quantity
getitem(T const & q, ssize_t i)
{
    return getitem(q, std::vector<ssize_t>{i});
}

template<typename T>
sycomore::Quantity const &
setitem(T & l, std::vector<ssize_t> const & i, sycomore::Quantity const & r)
{
    l[normalize_index(l.magnitude, i)] = r;
    return r;
}

template<typename T>
sycomore::Quantity const &
setitem(T & l, ssize_t i, sycomore::Quantity const & r)
{
    return setitem(l, std::vector<ssize_t>{i}, r);
}

#undef OUT_OF_PLACE_R_OPERATOR
#undef OUT_OF_PLACE_OPERATOR
#undef IN_PLACE_OPERATOR

}

}

#endif // _e0c9f20a_f96b_4c09_8459_f2be8e7df213
