#ifndef _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f
#define _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f

#include <ostream>
#include <xtensor/xtensor.hpp>

#include "sycomore/Array.h"
#include "sycomore/Dimensions.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

/// @brief Quantity in the SI system.
class SYCOMORE_API Quantity
{
public:
    /// @brief Magnitude of the quantity in base units.
    double magnitude;
    
    /// @brief Dimensions of the quantity.
    Dimensions dimensions;
    
    /// @brief Create a quantity from a magnitude and dimensions. 
    Quantity(double magnitude={}, Dimensions const & dimensions={});
    
    /// @brief Test whether magnitudes and dimensions are equal.
    bool operator==(Quantity const & other) const;
    
    /// @brief Test whether magnitudes or dimensions differ.
    bool operator!=(Quantity const & other) const;
    
    /// @brief In-place addition of a compatible quantity.
    Quantity & operator+=(Quantity const & other);
    
    /// @brief In-place addition of a compatible quantity.
    Quantity & operator+=(double s);
    
    /// @brief In-place subtraction of a compatible quantity.
    Quantity & operator-=(Quantity const & other);
    
    /// @brief In-place subtraction of a compatible quantity.
    Quantity & operator-=(double s);
    
    /// @brief In-place multiplication.
    Quantity & operator*=(Quantity const & other);
    
    /// @brief In-place multiplication.
    Quantity & operator*=(double scalar);
    
    /// @brief In-place division.
    Quantity & operator/=(Quantity const & other);
    
    /// @brief In-place division.
    Quantity & operator/=(double scalar);
    
    /// @brief In-place floating-point modulo.
    Quantity & operator%=(Quantity const & other);
    
    /// @brief In-place floating-point modulo.
    Quantity & operator%=(double scalar);

    /**
     * @brief Return the scalar value of the quantity converted to the given 
     * unit.
     *
     * Raise an exception if the given unit is not compatible.
    */
    double convert_to(Quantity const & destination) const;
    
    /**
     * @brief Convert to a scalar.
     *
     * Raise an exception if the quantity is not unitless.
     */
    operator double() const;
};

/// @brief Test whether magnitudes and dimensions are equal.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator==(Quantity const & q, T s) { return q == Quantity(s); }

/// @brief Test whether magnitudes and dimensions are equal.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator==(T s, Quantity const & q) { return q == Quantity(s); }

/// @brief Test whether magnitudes or dimensions differ.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator!=(Quantity const & q, T s) { return q != Quantity(s); }

/// @brief Test whether magnitudes or dimensions differ.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator!=(T s, Quantity const & q) { return q != Quantity(s); }

/// @brief Compare the magnitude of two compatible quantities.
SYCOMORE_API bool operator<(Quantity const & l, Quantity const & r);

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<(Quantity const & q, T s) { return q < Quantity(s); }

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<(T s, Quantity const & q) { return Quantity(s) < q; }

/// @brief Compare the magnitude of two compatible quantities.
SYCOMORE_API bool operator<=(Quantity const & l, Quantity const & r);

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<=(Quantity const & q, T s) { return q <= Quantity(s); }

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<=(T s, Quantity const & q) { return Quantity(s) <= q; }

/// @brief Compare the magnitude of two compatible quantities.
SYCOMORE_API bool operator>(Quantity const & l, Quantity const & r);

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>(Quantity const & q, T s) { return q > Quantity(s); }

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>(T s, Quantity const & q) { return Quantity(s) > q; }

/// @brief Compare the magnitude of two compatible quantities.
SYCOMORE_API bool operator>=(Quantity const & l, Quantity const & r);

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>=(Quantity const & q, T s) { return q >= Quantity(s); }

/// @brief Compare the magnitude of two compatible quantities.
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>=(T s, Quantity const & q) { return Quantity(s) >= q; }

/// @brief Identity operator.
SYCOMORE_API Quantity operator+(Quantity q);

/// @brief Return a quantity with the opposite magnitude
SYCOMORE_API Quantity operator-(Quantity q);

/// @brief Addition of compatible quantities
SYCOMORE_API Quantity operator+(Quantity l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>

/// @brief Addition of compatible quantities
Quantity operator+(Quantity const & q, T s) { return q+Quantity(s); }

/// @brief Addition of compatible quantities
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator+(T s, Quantity const & q) { return Quantity(s)+q; }

/// @brief Subtraction of compatible quantities
SYCOMORE_API Quantity operator-(Quantity l, Quantity const & r);

/// @brief Subtraction of compatible quantities
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator-(Quantity const & q, T s) { return q-Quantity(s); }

/// @brief Subtraction of compatible quantities
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator-(T s, Quantity const & q) { return Quantity(s)-q; }

/// @brief Multiplication
SYCOMORE_API Quantity operator*(Quantity l, Quantity const & r);

/// @brief Multiplication
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator*(Quantity q, T s) { q *= s; return q; }

/// @brief Multiplication
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator*(T s, Quantity const & q) { return q*s; }

/// @brief Division
SYCOMORE_API Quantity operator/(Quantity l, Quantity const & r);

/// @brief Division
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator/(Quantity q, T s) { q /= s; return q; }

/// @brief Division
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator/(T s, Quantity const & q) { return Quantity(s)/q; }

/// @brief Floating-point modulo
SYCOMORE_API Quantity operator%(Quantity l, Quantity const & r);

/// @brief Floating-point modulo
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator%(Quantity q, T s) { q %= double(s); return q; }

/// @brief String representation of a quantity
SYCOMORE_API std::ostream & operator<<(std::ostream & stream, Quantity const & q);

/// @brief Static-dimension array of Quantity objects
template<std::size_t N, xt::layout_type L=XTENSOR_DEFAULT_LAYOUT>
using TensorQ = xt::xtensor<Quantity, N, L>;

/// @brief Dynamic-dimension array of Quantity objects
using ArrayQ = xt::xarray<Quantity>;

/// @brief 3D vector of Quantity
using Vector3Q = Vector3<Quantity>;

/// @brief 3x3 matrix of Quantity
using Matrix3x3Q = Matrix3x3<Quantity>;

/// @brief Convert a sequence of Quantity to given unit
template<typename SourceIt, typename DestinationIt>
DestinationIt convert_to(
    SourceIt const & begin, SourceIt const & end, DestinationIt destination,
    Quantity const & t)
{
    return std::transform(
        begin, end, destination, [&](auto && x) { return x.convert_to(t); });
}

/// @brief Convert a sequence of Quantity to given unit
template<std::size_t N, xt::layout_type L=XTENSOR_DEFAULT_LAYOUT>
TensorR<N, L> convert_to(TensorQ<N, L> const & q, Quantity const & t)
{
    TensorR<N, L> r(q.shape());
    convert_to(q.begin(), q.end(), r.begin(), t);
    return r;
}

}

namespace std
{

/// @brief Return a quantity with the absolute value of the magnitude
SYCOMORE_API sycomore::Quantity abs(sycomore::Quantity q);

/// @brief Raise a quantity to a power
SYCOMORE_API sycomore::Quantity pow(sycomore::Quantity q, double e);

/// @brief Round the magnitude of a quantity
SYCOMORE_API sycomore::Quantity round(sycomore::Quantity q);

/// @brief Truncate the magnitude of a quantity
SYCOMORE_API sycomore::Quantity trunc(sycomore::Quantity q);

/// @brief Quantity with the largest integer magnitude not greater than the magnitude
SYCOMORE_API sycomore::Quantity floor(sycomore::Quantity q);

/// @brief Quantity with the smallest integer magnitude not less than the magnitude
SYCOMORE_API sycomore::Quantity ceil(sycomore::Quantity q);

/// @brief Hash functor
template<>
struct SYCOMORE_API hash<sycomore::Quantity>
{
    /// @brief Hash function
    std::size_t operator()(sycomore::Quantity const & q) const;
};

}

#endif // _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f
