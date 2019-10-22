#include "Array.h"

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <numeric>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace sycomore
{

template<typename T>
Array<T>
::Array()
: Array(0)
{
    // Nothing else.
}

template<typename T>
Array<T>
::Array(std::size_t count)
: _size(count), _is_view(false), _data(nullptr)
{
    this->_data = new T[this->_size];
}

template<typename T>
Array<T>
::Array(std::size_t count, T const & value)
: Array(count)
{
    std::fill(this->begin(), this->end(), value);
}

template<typename T>
Array<T>
::Array(T * pointer, std::size_t size)
: _size(size), _is_view(true), _data(pointer)
{
    // Nothing else.
}

template<typename T>
Array<T>
::Array(Array<T> const & other)
: Array(other._size)
{
    std::copy(other.begin(), other.end(), this->begin());
}

template<typename T>
Array<T>
::Array(std::initializer_list<T> const & initializer_list)
: Array(initializer_list.size())
{
    std::copy(initializer_list.begin(), initializer_list.end(), this->begin());
}

template<typename T>
Array<T>
::Array(Array<T> && other)
: _size(other._size), _is_view(other._is_view), _data(other._data)
{
    other._size = 0;
    other._is_view = false;
    other._data = nullptr;
}

template<typename T>
Array<T> &
Array<T>
::operator=(Array<T> const & other)
{
    if(this->_is_view)
    {
        if(this->_size != other.size())
        {
            throw std::runtime_error("Assigning to a view must preserve size");
        }
    }
    else
    {
        this->_size = other._size;

        if(this->_data != nullptr)
        {
            delete[] this->_data;
        }
        this->_data = new T[this->_size];
    }
    std::copy(other.begin(), other.end(), this->begin());
    return *this;
}

template<typename T>
Array<T> &
Array<T>
::operator=(Array<T> && other)
{
    this->_size = other._size;
    this->_is_view = other._is_view;
    this->_data = other._data;

    other._size = 0;
    other._is_view = false;
    other._data = nullptr;

    return *this;
}

template<typename T>
Array<T>
::~Array()
{
    if(!this->_is_view && this->_data != nullptr)
    {
        delete[] this->_data;
        this->_data = nullptr;
        this->_size = 0;
    }
}

template<typename T>
std::size_t
Array<T>
::size() const
{
    return this->_size;
}

template<typename T>
bool
Array<T>
::empty() const
{
    return this->size() == 0;
}

template<typename T>
bool
Array<T>
::is_view() const
{
    return this->_is_view;
}

template<typename T>
T const *
Array<T>
::data() const
{
    return this->_data;
}

template<typename T>
T *
Array<T>
::data()
{
    return this->_data;
}

template<typename T>
template<typename T2>
Array<T2>
Array<T>
::astype() const
{
    Array<T2> result(this->size());
    auto source = this->data();
    auto destination = result.data();
    for(std::size_t i=0; i<this->size(); ++i)
    {
        *destination = static_cast<T>(*source);
        ++source;
        ++destination;
    }
    return result;
}

template<typename T>
T const &
Array<T>
::operator[](std::size_t i) const
{
    return this->_data[i];
}

template<typename T>
T &
Array<T>
::operator[](std::size_t i)
{
    return this->_data[i];
}

#define SYCOMORE_MUTATING_OPERATOR_SCALAR(op) \
template<typename T> \
Array<T> & \
Array<T> \
::operator op##=(value_type const & s) \
{ \
    for(std::size_t i=0; i!=this->_size; ++i) \
    { \
        this->_data[i] op##= s; \
    } \
    return *this; \
}

SYCOMORE_MUTATING_OPERATOR_SCALAR(+)
SYCOMORE_MUTATING_OPERATOR_SCALAR(-)
SYCOMORE_MUTATING_OPERATOR_SCALAR(*)
SYCOMORE_MUTATING_OPERATOR_SCALAR(/)

#undef SYCOMORE_MUTATING_OPERATOR_SCALAR

#define SYCOMORE_MUTATING_OPERATOR_ARRAY(op) \
template<typename T> \
Array<T> & \
Array<T> \
::operator op##=(Array<T> const & other) \
{ \
    if(this->size() != other.size()) \
    { \
        throw std::runtime_error("Size mismatch"); \
    } \
    for(std::size_t i=0; i!=this->_size; ++i) \
    { \
        this->_data[i] op##= other._data[i]; \
    } \
    return *this; \
}

SYCOMORE_MUTATING_OPERATOR_ARRAY(+)
SYCOMORE_MUTATING_OPERATOR_ARRAY(-)

#undef SYCOMORE_MUTATING_OPERATOR_ARRAY

template<typename T>
typename Array<T>::const_iterator
Array<T>
::begin() const
{
    return this->_data;
}

template<typename T>
typename Array<T>::const_iterator
Array<T>
::cbegin() const
{
    return this->_data;
}

template<typename T>
typename Array<T>::iterator
Array<T>
::begin()
{
    return this->_data;
}

template<typename T>
typename Array<T>::const_iterator
Array<T>
::end() const
{
    return this->_data+this->_size;
}

template<typename T>
typename Array<T>::const_iterator
Array<T>
::cend() const
{
    return this->_data+this->_size;
}

template<typename T>
typename Array<T>::iterator
Array<T>
::end()
{
    return this->_data+this->_size;
}

template<typename T>
bool operator==(Array<T> const & l, Array<T> const & r)
{
    bool equal=true;
    if(l.size() != r.size())
    {
        equal = false;
    }
    else
    {
        for(std::size_t i=0; i!=l.size() && equal; ++i)
        {
            if(l[i] != r[i])
            {
                equal = false;
            }
        }
    }
    return equal;
}

template<typename T>
bool operator!=(Array<T> const & l, Array<T> const & r)
{
    return !(l==r);
}

template<typename T>
bool operator<(Array<T> const & l, Array<T> const & r)
{
    return std::lexicographical_compare(l.begin(), l.end(), r.begin(), r.end());
}

template<typename T>
bool operator>(Array<T> const & l, Array<T> const & r)
{
    return std::lexicographical_compare(
        l.begin(), l.end(), r.begin(), r.end(), std::greater<T>());
}

template<typename T>
bool operator<=(Array<T> const & l, Array<T> const & r)
{
    return !(l > r);
}

template<typename T>
bool operator>=(Array<T> const & l, Array<T> const & r)
{
    return !(l < r);
}

template<typename T>
Array<T> operator-(Array<T> l)
{
    for(std::size_t i=0; i!=l.size(); ++i)
    {
        l[i] = -l[i];
    }
    return l;
}

#define SYCOMORE_OPERATOR_SCALAR_RIGHT(op) \
template<typename T1, typename T2> \
Array<decltype(std::declval<T1>() op std::declval<T2>())> \
operator op(Array<T1> const & a, T2 const & s) \
{ \
    Array<decltype(std::declval<T1>() op std::declval<T2>())> result(a.size()); \
    for(std::size_t i=0; i<result.size(); ++i) \
    { \
        result[i] = a[i] op s; \
    } \
    return result; \
}

SYCOMORE_OPERATOR_SCALAR_RIGHT(+);
SYCOMORE_OPERATOR_SCALAR_RIGHT(-);
SYCOMORE_OPERATOR_SCALAR_RIGHT(*);
SYCOMORE_OPERATOR_SCALAR_RIGHT(/);

#undef SYCOMORE_OPERATOR_SCALAR_RIGHT

#define SYCOMORE_OPERATOR_SCALAR_LEFT(op) \
template<typename T1, typename T2> \
Array<decltype(std::declval<T1>() op std::declval<T2>())> \
operator op(T1 const & s, Array<T2> const & a) \
{ \
    Array<decltype(std::declval<T1>() op std::declval<T2>())> result(a.size()); \
    for(std::size_t i=0; i<result.size(); ++i) \
    { \
        result[i] = a[i] op s; \
    } \
    return result; \
}

SYCOMORE_OPERATOR_SCALAR_LEFT(+);
SYCOMORE_OPERATOR_SCALAR_LEFT(*);

#undef SYCOMORE_OPERATOR_SCALAR_LEFT

#define SYCOMORE_OPERATOR_ARRAY(op) \
template<typename T1, typename T2> \
Array<decltype(std::declval<T1>() op std::declval<T2>())> \
operator op(Array<T1> const & l, Array<T2> const & r) \
{ \
    Array<decltype(std::declval<T1>() op std::declval<T2>())> x(std::min(l.size(), r.size())); \
    for(std::size_t i=0; i!=x.size(); ++i) \
    { \
        x[i] = l[i] op r[i]; \
    } \
    return x; \
}

SYCOMORE_OPERATOR_ARRAY(+)
SYCOMORE_OPERATOR_ARRAY(-)

#undef SYCOMORE_OPERATOR_ARRAY

template<typename T>
Array<T> minimum(Array<T> const & l, Array<T> const & r)
{
    Array<T> result(std::min(l.size(), r.size()));
    for(std::size_t i=0; i!=result.size(); ++i)
    {
        result[i] = std::min(l[i], r[i]);
    }
    return result;
}

template<typename T>
Array<T> maximum(Array<T> const & l, Array<T> const & r)
{
    Array<T> result(std::min(l.size(), r.size()));
    for(std::size_t i=0; i!=result.size(); ++i)
    {
        result[i] = std::max(l[i], r[i]);
    }
    return result;
}

template<typename T1, typename T2>
decltype(std::declval<T1>()*std::declval<T2>())
dot(Array<T1> const & l, Array<T2> const & r)
{
    decltype(std::declval<T1>()*std::declval<T2>()) result(0);
    auto const size=std::min(l.size(), r.size());
    for(std::size_t i=0; i!=size; ++i)
    {
        result += l[i]*r[i];
    }
    return result;
}

template<typename T>
std::ostream & operator<<(std::ostream & stream, Array<T> const & array)
{
    stream << "(";
    if(!array.empty())
    {
        auto printer = [](T const & x) {
            std::stringstream s; s << x; return s.str(); };
        auto const string = std::accumulate(
            std::next(array.begin()), array.end(), printer(*array.begin()),
            [&](std::string const & s, T const & x) { return s+" "+printer(x); });
        stream << string;
    }
    stream << ")";
    return stream;
}

}

namespace std
{

template<typename T>
sycomore::Array<T> abs(sycomore::Array<T> const & a)
{
    sycomore::Array<T> result(a);
    for(auto & x: result)
    {
        x = std::abs(x);
    }
    return result;
}

template<typename T>
std::size_t
hash<sycomore::Array<T>>
::operator()(sycomore::Array<T> const & index) const
{
    return sycomore::hash_range(index.cbegin(), index.cend());
}

}
