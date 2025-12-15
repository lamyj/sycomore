#ifndef _d81aa20b_dd7e_4a5a_a47c_837fd494e686
#define _d81aa20b_dd7e_4a5a_a47c_837fd494e686

#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include "sycomore/QuantityArray.h"
#include "sycomore/Quantity.h"
#include "sycomore/QuantityConstView.h"

// NOTE: forbid TensorQ to participate in the operators overloads of namespace xt
namespace xt
{
namespace detail
{
template<typename F, typename T>
struct xfunction_type<F, T, decltype(pybind11::self)> { };
template<typename F, typename T>
struct xfunction_type<F, decltype(pybind11::self), T> { };
template<typename F, typename T>
struct xfunction_type<F, T, decltype(pybind11::self) const &> { };
template<typename F, typename T>
struct xfunction_type<F, decltype(pybind11::self) const &, T> { };
}
}

namespace sycomore
{

namespace wrappers
{

template<typename T>
pybind11::class_<T>
wrap_quantity_class(pybind11::module & m, std::string const & name);

template<typename T>
pybind11::class_<T>
wrap_quantity_array(pybind11::module & m, pybind11::class_<T> & _class);

template<typename T>
T as_quantity(pybind11::array_t<pybind11::object> array);

std::size_t normalize_index(std::size_t shape, ssize_t i);

template<typename T>
std::vector<std::size_t>
normalize_index(T const & magnitude, std::vector<ssize_t> const & i);

template<typename T>
pybind11::object getitem(T const & q, pybind11::object index);

template<typename T>
sycomore::Quantity const &
setitem(T & l, std::vector<ssize_t> const & i, sycomore::Quantity const & r);

template<typename T>
sycomore::Quantity const &
setitem(T & l, ssize_t i, sycomore::Quantity const & r);

template<typename T>
struct QuantityConstIteratorAdapter
{
    /// @brief C++ object to iterate on
    T const & q;
    /// @brief Keep-alive reference to Python object
    pybind11::object r;
    /// @brief Index of current item
    std::size_t index;
    
    QuantityConstIteratorAdapter(T const & q, pybind11::object r);
    
    pybind11::object next();
};

}

}

#include "Quantity.txx"

#endif // _d81aa20b_dd7e_4a5a_a47c_837fd494e686
