#ifndef _950e1518_e22c_46bd_afd1_5f2d045d8d73
#define _950e1518_e22c_46bd_afd1_5f2d045d8d73

#include <cmath>
#include <stdexcept>

#include <xtensor/xmath.hpp>

#include "sycomore/Dimensions.h"
#include "sycomore/quantity_traits.h"

namespace sycomore
{

/**
 * @brief Common semantics class for all quantity-like objects
 *
 * This follows the CRTP pattern. The magnitude type of derived must be either
 * a built-in floating point type or and xtensor container.
 */
template<typename TDerived>
class QuantityInterface
{
public:
    QuantityInterface() = default;
    QuantityInterface(QuantityInterface<TDerived> const &) = default;
    QuantityInterface(QuantityInterface<TDerived> &&) = default;
    QuantityInterface<TDerived> & operator=(
        QuantityInterface<TDerived> const &) = default;
    QuantityInterface<TDerived> & operator=(
        QuantityInterface<TDerived> &&) = default;
    virtual ~QuantityInterface() = default;
    
    /**
     * @brief Raise an exception if this->dimensions do not match provided
     * dimensions
     * @param message message associated with the exception
     */
    void check_dimensions(Dimensions const & dimensions, std::string const & message="") const;

    /**
     * @brief Raise an exception if this->dimensions do not match the dimensions
     * of the provided object
     * @param message message associated with the exception
     */
    template<typename T>
    void check_dimensions(T const & other, std::string const & message="") const;
    
    /// @brief Return the concrete object
    virtual TDerived & derived_cast();
    
    /// @brief Return the concrete object
    virtual TDerived const & derived_cast() const;
    
    /// @brief Test whether magnitudes and dimensions are equal
    bool operator==(TDerived const & right) const;
    
    /// @brief Test whether magnitudes or dimensions differ
    bool operator!=(TDerived const & right) const;
    
    /// @brief In-place addition of a compatible quantity
    template<typename T>
    TDerived & operator+=(T const & right);
    
    /// @brief In-place subtraction of a compatible quantity
    template<typename T>
    TDerived & operator-=(T const & right);
    
    /// @brief In-place multiplication
    template<typename T, enable_if_quantity<T> = true>
    TDerived & operator*=(T const & right)
    {
        auto & left = this->derived_cast();
        left.magnitude *= right.magnitude;
        left.dimensions *= right.dimensions;
        return left;
    }
    
    /// @brief In-place multiplication
    template<typename T, enable_if_not_quantity<T> = true>
    TDerived & operator*=(T const & right)
    {
        auto & left = this->derived_cast();
        left.magnitude *= right;
        return left;
    }
    
    /// @brief In-place division
    template<typename T, enable_if_quantity<T> = true>
    TDerived & operator/=(T const & right)
    {
        auto & left = this->derived_cast();
        left.magnitude /= right.magnitude;
        left.dimensions /= right.dimensions;
        return left;
    }
    
    /// @brief In-place division
    template<typename T, enable_if_not_quantity<T> = true>
    TDerived & operator/=(T const & right)
    {
        auto & left = this->derived_cast();
        left.magnitude /= right;
        return left;
    }
    
    /// @brief In-place floating-point modulo
    template<typename T, enable_if_quantity<T> = true>
    TDerived & operator%=(T const & right)
    {
        auto & left = this->derived_cast();
        if(right.dimensions == Dimensionless || left.dimensions == right.dimensions)
        {
            return this->operator%=(right.magnitude);
        }
        else
        {
            throw std::runtime_error("Modulo requires same dimensions");
        }
        return left;
    }
    
    /// @brief In-place floating-point modulo
    template<typename T, enable_if_not_quantity<T> = true>
    TDerived & operator%=(T const & right)
    {
        auto & left = this->derived_cast();
        using std::fmod;
        using xt::fmod;
        left.magnitude = fmod(left.magnitude, right);
        return left;
    }
    
    template<typename T>
    TDerived & fill(T const & x);
};

/// @brief Explicit specialization for true quantity types
template<typename T>
struct is_quantity<
    T,
    std::enable_if_t<std::is_base_of<QuantityInterface<T>, T>::value>>
: public std::true_type { };

/// @brief Identity operator
template<typename T, enable_if_quantity<T> = true>
OwningType<T> operator+(T x)
{
    return {+x.magnitude, x.dimensions};
}

/// @brief Return a quantity with the opposite magnitude
template<typename T, enable_if_quantity<T> = true>
OwningType<T> operator-(T const & x)
{
    return {-x.magnitude, x.dimensions};
}

/// @brief Addition of compatible quantities
template<typename T1, typename T2, enable_if_quantity<T1> = true, enable_if_quantity<T2> = true>
CommonQuantityType<T1, T2> operator+(T1 const & left, T2 const & right)
{
    left.check_dimensions(right);
    return {left.magnitude + right.magnitude, left.dimensions};
}

/// @brief Subtraction of compatible quantities
template<typename T1, typename T2, enable_if_quantity<T1> = true, enable_if_quantity<T2> = true>
CommonQuantityType<T1, T2> operator-(T1 const & left, T2 const & right)
{
    left.check_dimensions(right);
    return {left.magnitude - right.magnitude, left.dimensions};
}

/// @brief Multiplication of quantities
template<typename T1, typename T2, enable_if_quantity<T1> = true, enable_if_quantity<T2> = true>
CommonQuantityType<T1, T2> operator*(T1 const & left, T2 const & right)
{
    return {left.magnitude * right.magnitude, left.dimensions * right.dimensions};
}

/// @brief Multiplication of a quantity and a scalar
template<typename T1, typename T2, enable_if_quantity<T1> = true, enable_if_not_quantity<T2> = true>
OwningType<T1> operator*(T1 const & left, T2 const & right)
{
    return {left.magnitude * right, left.dimensions};
}

/// @brief Multiplication of a scalar and a quantity
template<typename T1, typename T2, enable_if_not_quantity<T1> = true, enable_if_quantity<T2> = true>
OwningType<T2> operator*(T1 const & left, T2 const & right)
{
    return {left * right.magnitude, right.dimensions};
}

/// @brief Division of quantities
template<typename T1, typename T2, enable_if_quantity<T1> = true, enable_if_quantity<T2> = true>
CommonQuantityType<T1, T2> operator/(T1 const & left, T2 right)
{
    return {left.magnitude / right.magnitude, left.dimensions / right.dimensions};
}

/// @brief Division of a quantity and a scalar
template<typename T1, typename T2, enable_if_quantity<T1> = true, enable_if_not_quantity<T2> = true>
OwningType<T1> operator/(T1 left, T2 const & right)
{
    return {left.magnitude / right, left.dimensions};
}

/// @brief Division of quantities
template<typename T1, typename T2, enable_if_not_quantity<T1> = true, enable_if_quantity<T2> = true>
OwningType<T2> operator/(T1 const & left, T2 right)
{
    return {left/right.magnitude, std::pow(right.dimensions, -1)};
}

/// @brief Floating point modulo of quantities
template<typename T1, typename T2, enable_if_quantity<T1> = true>
OwningType<T1> operator%(T1 left, T2 const & right)
{
    return OwningType<T1>{left.magnitude, left.dimensions} %= right;
}

}

namespace std
{

/// @brief Return a quantity with the absolute value of the magnitude
template<typename T, sycomore::enable_if_quantity<T> = true>
sycomore::OwningType<T> abs(T const & q)
{
    using std::abs;
    using xt::abs;
    return {abs(q.magnitude), q.dimensions};
}

/// @brief Raise a quantity to a power
template<typename T, sycomore::enable_if_quantity<T> = true>
sycomore::OwningType<T> pow(T const & q, double e)
{
    using std::pow;
    using xt::pow;
    return {pow(q.magnitude, e), pow(q.dimensions, e)};
}

/// @brief Round the magnitude of a quantity
template<typename T, sycomore::enable_if_quantity<T> = true>
sycomore::OwningType<T> round(T const & q)
{
    using std::round;
    using xt::round;
    return {round(q.magnitude), q.dimensions};
}

/// @brief Truncate the magnitude of a quantity
template<typename T, sycomore::enable_if_quantity<T> = true>
sycomore::OwningType<T> trunc(T const & q)
{
    using std::trunc;
    using xt::trunc;
    return {trunc(q.magnitude), q.dimensions};
}

/// @brief Quantity with the largest integer magnitude not greater than the magnitude
template<typename T, sycomore::enable_if_quantity<T> = true>
sycomore::OwningType<T> floor(T const & q)
{
    using std::floor;
    using xt::floor;
    return {floor(q.magnitude), q.dimensions};
}

/// @brief Quantity with the smallest integer magnitude not less than the magnitude
template<typename T, sycomore::enable_if_quantity<T> = true>
sycomore::OwningType<T> ceil(T const & q)
{
    using std::ceil;
    using xt::ceil;
    return {ceil(q.magnitude), q.dimensions};
}

}

#include "QuantityInterface.txx"

#endif // _950e1518_e22c_46bd_afd1_5f2d045d8d73
