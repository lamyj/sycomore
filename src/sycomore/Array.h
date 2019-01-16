#ifndef _dc7a30fd_6048_4dfc_a551_efae434269f0
#define _dc7a30fd_6048_4dfc_a551_efae434269f0

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <type_traits>
#include <utility>

#include <iostream>

namespace sycomore
{

/// @brief One-dimensional array with arithmetic operations.
template<typename T>
class Array
{
public:
    using value_type = T;

    explicit Array(size_t count=0);
    Array(size_t count, T const & value);
    Array(Array<T> const & other);
    Array(std::initializer_list<T> const & init);

    Array(Array<T> && other);

    Array<T> & operator=(Array<T> const & other);

    Array<T> & operator=(Array<T> && other);

    ~Array();

    /// @brief Return the number of elements.
    size_t size() const;

    /// @brief Test whether the array contains no element.
    bool empty() const;

    /// @brief Read-only access, no bounds checking.
    T const & operator[](size_t i) const;

    /// @brief Read-write access, no bounds checking.
    T & operator[](size_t i);

    /// @addtogroup mutating_operators Mutating operators
    /// @{

    Array<T> & operator+=(value_type const & s);

    Array<T> & operator-=(value_type const & s);

    Array<T> & operator*=(value_type const & s);

    Array<T> & operator/=(value_type const & s);

    Array<T> & operator+=(Array<T> const & other);

    Array<T> & operator-=(Array<T> const & other);

    /// @}

    /// @addtogroup iterators Iterators
    /// @{
    using iterator = T *;
    using const_iterator = T const *;

    const_iterator begin() const;
    const_iterator cbegin() const;
    iterator begin();

    const_iterator end() const;
    const_iterator cend() const;
    iterator end();
    /// @}
private:
    /// @brief Number of elements in the array
    size_t _size;

    /// @brief Contiguous memory area holding the elements
    T * _data;
};

template<typename T>
bool operator==(Array<T> const & l, Array<T> const & r);

template<typename T>
bool operator!=(Array<T> const & l, Array<T> const & r);

template<typename T>
bool operator<(Array<T> const & l, Array<T> const & r);

template<typename T>
bool operator>(Array<T> const & l, Array<T> const & r);

template<typename T>
bool operator<=(Array<T> const & l, Array<T> const & r);

template<typename T>
bool operator>=(Array<T> const & l, Array<T> const & r);

template<typename T>
Array<T> operator-(Array<T> l);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()+std::declval<T2>())>
operator+(Array<T1> const & l, Array<T2> const & r);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()-std::declval<T2>())>
operator-(Array<T1> const & l, Array<T2> const & r);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()+std::declval<T2>())>
operator+(Array<T1> const & l, T2 const & s);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()+std::declval<T2>())>
operator+(T1 const & s, Array<T2> const & l);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()-std::declval<T2>())>
operator-(Array<T1> const & l, T2 const & s);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()*std::declval<T2>())>
operator*(Array<T1> const & l, T2 const & s);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()*std::declval<T2>())>
operator*(T1 const & s, Array<T2> const &l);

template<typename T1, typename T2>
Array<decltype(std::declval<T1>()/std::declval<T2>())>
operator/(Array<T1> const & l, T2 const & s);

template<typename T>
Array<T> minimum(Array<T> const & l, Array<T> const & r);

template<typename T>
Array<T> maximum(Array<T> const & l, Array<T> const & r);

template<typename T1, typename T2>
decltype(std::declval<T1>()*std::declval<T2>())
dot(Array<T1> const & l, Array<T2> const & r);

template<typename T>
std::ostream & operator<<(std::ostream & stream, Array<T> const & index);

}

#include "Array.txx"

#endif // _dc7a30fd_6048_4dfc_a551_efae434269f0
