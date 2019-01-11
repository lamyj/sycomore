#include "Array.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <type_traits>
#include <utility>

namespace sycomore
{

template<typename T>
Array<T>
::Array(size_t count)
: _size(count), _data(nullptr)
{
    this->_data = new T[this->_size];
}

template<typename T>
Array<T>
::Array(size_t count, T const & value)
: Array(count)
{
    std::fill(this->begin(), this->end(), value);
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
: _size(other._size), _data(other._data)
{
    other._size = 0;
    other._data = nullptr;
}

//template<typename T>
//template<typename Expression, typename std::enable_if<ArrayExpressionTraits<Expression>::is_array_expression, int>::type>
//Array<T>
//::Array(Expression const & expression)
//: Array(expression.size())
//{
//    for(size_t i=0; i<expression.size(); ++i)
//    {
//        this->_data[i] = expression[i];
//    }
//}

template<typename T>
Array<T> &
Array<T>
::operator=(Array<T> const & other)
{
    this->_size = other._size;

    if(this->_data != nullptr)
    {
        delete[] this->_data;
    }
    this->_data = new T[this->_size];

    std::copy(other.begin(), other.end(), this->begin());
    return *this;
}

template<typename T>
Array<T> &
Array<T>
::operator=(Array<T> && other)
{
    this->_size = other._size;
    this->_data = other._data;

    other._size = 0;
    other._data = nullptr;

    return *this;
}

//template<typename T>
//template<typename Expression, typename std::enable_if<ArrayExpressionTraits<Expression>::is_array_expression, int>::type>
//Array<T> &
//Array<T>
//::operator=(Expression const & expression)
//{
//    this->_size = expression.size();

//    if(this->_data != nullptr)
//    {
//        delete[] this->_data;
//    }
//    this->_data = new T[this->_size];

//    for(size_t i=0; i<expression.size(); ++i)
//    {
//        this->_data[i] = expression[i];
//    }

//    return *this;
//}

template<typename T>
Array<T>
::~Array()
{
    if(this->_data != nullptr)
    {
        delete[] this->_data;
        this->_data = nullptr;
        this->_size = 0;
    }
}

template<typename T>
size_t
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
T const &
Array<T>
::operator[](size_t i) const
{
    return this->_data[i];
}

template<typename T>
T &
Array<T>
::operator[](size_t i)
{
    return this->_data[i];
}

#define SYCOMORE_MUTATING_OPERATOR_SCALAR(op) \
template<typename T> \
Array<T> & \
Array<T> \
::operator op##=(value_type const & s) \
{ \
    for(size_t i=0; i!=this->_size; ++i) \
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
    for(size_t i=0; i!=this->_size; ++i) \
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
        for(size_t i=0; i!=l.size() && equal; ++i)
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
    for(size_t i=0; i!=l.size(); ++i)
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
    Array<decltype(std::declval<T1>() op std::declval<T2>())> result; \
    for(size_t i=0; i<result.size(); ++i) \
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
    Array<decltype(std::declval<T1>() op std::declval<T2>())> result; \
    for(size_t i=0; i<result.size(); ++i) \
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
    for(size_t i=0; i!=x.size(); ++i) \
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
    for(size_t i=0; i!=result.size(); ++i)
    {
        result[i] = std::min(l[i], r[i]);
    }
    return result;
}

template<typename T>
Array<T> maximum(Array<T> const & l, Array<T> const & r)
{
    Array<T> result(std::min(l.size(), r.size()));
    for(size_t i=0; i!=result.size(); ++i)
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
    for(size_t i=0; i!=size; ++i)
    {
        result += l[i]*r[i];
    }
    return result;
}

}
