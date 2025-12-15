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
        .def_property_readonly("scalar", &T::scalar)
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
    /********************* (OUT OF PLACE, HETEROGENENOUS) *********************/
    /**************************************************************************/
    _class
        .def(self + Container()).def(Container() + self)
        .def(self - Container()).def(Container() - self)
        .def(self * Container()).def(Container() * self)
        .def(self / Container()).def(Container() / self)
        .def(
            "__floordiv__",
            [](T const & l, Container const & r) { return std::floor(l/r); },
            is_operator())
        .def(
            "__rfloordiv__",
            [](T const & r, Container const & l) { return std::floor(l/r); },
            is_operator())
        .def("__mod__", sycomore::fmod<T, Container>, is_operator())
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
                "magnitude");
    
    return _class;
}

template<typename T>
pybind11::class_<T>
wrap_quantity_array(pybind11::module & m, pybind11::class_<T> & _class)
{
    using namespace pybind11;
    using namespace sycomore;
    
    _class
        .def(init(&wrappers::as_quantity<T>))
        .def(self += Quantity()).def(self + Quantity()).def(Quantity() + self)
        .def(self -= Quantity()).def(self - Quantity()).def(Quantity() - self)
        .def(self *= Quantity()).def(self * Quantity()).def(Quantity() * self)
        .def(self /= Quantity()).def(self / Quantity()).def(Quantity() / self)
        .def(
            "__ifloordiv__",
            [](T & l, Quantity const & r) { l = std::floor(l/r); return l; },
            is_operator())
        .def(
            "__floordiv__",
            [](T const & l, Quantity const & r) { return std::floor(l/r); },
            is_operator())
        .def(
            "__rfloordiv__",
            [](T const & r, Quantity const & l) { return std::floor(l/r); },
            is_operator())
        .def(
            "__imod__",
            [](T & l, Quantity const & r) { return (l %= r); },
            is_operator());
    
    pybind11::class_<QuantityConstIteratorAdapter<T>>(
            m, (std::string(typeid(T).name())+"Iterator").c_str())
        .def("__iter__", [](QuantityConstIteratorAdapter<T> it) { return it; })
        .def("__next__", &QuantityConstIteratorAdapter<T>::next);
    
    _class
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
        .def("__getitem__", &getitem<T>)
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
            "__iter__", 
            [](object q) {
                return QuantityConstIteratorAdapter<T>(q.cast<T const &>(), q);
            }
        );
    
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
            "Number of arguments ("
            + std::to_string(i.size())
            + ") does not match the number of dimensions ("
            + std::to_string(magnitude.dimension())
            + ")");
    }
    
    // Convert signed values (counting from end) to unsigned values
    std::vector<std::size_t> unsigned_i(i.size());
    for(std::size_t d=0; d != i.size(); ++d)
    {
        if(i[d] < 0)
        {
            unsigned_i[d] = magnitude.shape()[d] + i[d];
        }
        else
        {
            unsigned_i[d] = i[d];
        }
    }
    
    xt::check_element_index(magnitude.shape(), unsigned_i.begin(), unsigned_i.end());
    
    return unsigned_i;
}

template<typename T>
pybind11::object
getitem(T const & q, pybind11::object index)
{
    if(pybind11::isinstance<pybind11::int_>(index))
    {
        auto const view = sycomore::view(
            q, normalize_index(q.shape(0), index.cast<ssize_t>()));
        return view.size() == 1
            ? pybind11::cast(sycomore::Quantity(view.unchecked(0)))
            : pybind11::cast(sycomore::ArrayQ(view));
    }
    else 
    {
        auto const length = pybind11::len(index);
        
        if(length > q.shape().size())
        {
            throw std::out_of_range(
                "Number of arguments ("
                + std::to_string(length)
                + ") is greater than the number of dimensions ("
                + std::to_string(q.shape().size()) + ")");
        }
        
        xt::xstrided_slice_vector slices;
        slices.reserve(length);
        
        // Dimension counter
        std::size_t d = 0;
        for(auto && item: index)
        {
            auto const s = q.shape(d);
            
            if(pybind11::isinstance<pybind11::int_>(item))
            {
                slices.push_back(normalize_index(s, item.cast<ssize_t>()));
            }
            else
            {
                auto const slice = item.cast<pybind11::slice>();
                auto const start = slice.attr("start").cast<pybind11::object>();
                auto const stop = slice.attr("stop").cast<pybind11::object>();
                auto const step = slice.attr("step").cast<pybind11::object>();
                
                auto const start_is_none = start.is(pybind11::none());
                auto const stop_is_none = stop.is(pybind11::none());
                auto const step_is_none = step.is(pybind11::none());
                
                // FIXME: is there a cleaner way to write all 8 different calls?
                // Note that xt::range returns a type which depends on the
                // parameters
                using namespace xt::placeholders;
                if(start_is_none && stop_is_none && step_is_none)
                {
                    slices.push_back(xt::all());
                }
                else if(start_is_none && stop_is_none && !step_is_none)
                {
                    slices.push_back(xt::range(_, _, step.cast<ssize_t>()));
                }
                else if(start_is_none && !stop_is_none && step_is_none)
                {
                    slices.push_back(xt::range(_, stop.cast<ssize_t>(), _));
                }
                else if(start_is_none && !stop_is_none && !step_is_none)
                {
                    slices.push_back(xt::range(_, stop.cast<ssize_t>(), step.cast<ssize_t>()));
                }
                else if(!start_is_none && stop_is_none && step_is_none)
                {
                    slices.push_back(xt::range(start.cast<ssize_t>(), _, _));
                }
                else if(!start_is_none && stop_is_none && !step_is_none)
                {
                    slices.push_back(xt::range(start.cast<ssize_t>(), _, step.cast<ssize_t>()));
                }
                else if(!start_is_none && !stop_is_none && step_is_none)
                {
                    slices.push_back(xt::range(start.cast<ssize_t>(), stop.cast<ssize_t>(), _));
                }
                else if(!start_is_none && !stop_is_none && !step_is_none)
                {
                    slices.push_back(xt::range(start.cast<ssize_t>(), stop.cast<ssize_t>(), step.cast<ssize_t>()));
                }
            }
            
            ++d;
        }
        
        auto const view = sycomore::strided_view(q, slices);
        return view.size() == 1
            ? pybind11::cast(sycomore::Quantity(view.unchecked(0)))
            : pybind11::cast(sycomore::ArrayQ(view));
    }
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

template<typename T>
QuantityConstIteratorAdapter<T>
::QuantityConstIteratorAdapter(T const & q, pybind11::object r)
: q(q), r(r), index(0)
{
    // Nothing else
}

template<typename T>
pybind11::object
QuantityConstIteratorAdapter<T>
::next()
{
    auto const & shape = this->q.shape();
    
    if(shape.empty() || this->index == shape[0])
    {
        throw pybind11::stop_iteration();
    }
    else
    {
        if(shape.size() == 1)
        {
            return pybind11::cast(sycomore::Quantity(this->q.unchecked(this->index++)));
        }
        else
        {
            return pybind11::cast(
                sycomore::ArrayQ(sycomore::view(this->q, this->index++)));
        }
    }
}

}

}

#endif // _e0c9f20a_f96b_4c09_8459_f2be8e7df213
