#ifndef _dfbc0517_611a_4989_a51c_fa60b94c587f
#define _dfbc0517_611a_4989_a51c_fa60b94c587f

#include <xtensor/xtensor.hpp>

#include "sycomore/hash.h"
#include "sycomore/QuantityBase.h"

namespace sycomore
{

/// @brief Scalar quantity
class Quantity: public QuantityBase<Quantity, double>
{
public:
    using Self = Quantity;
    using Container = double;
    using Base = QuantityBase<Self, Container>;
    
    using Base::Base;
    
    Quantity(Quantity const &) = default;
    Quantity(Quantity &&) = default;
    Quantity & operator=(Quantity const &) = default;
    Quantity & operator=(Quantity &&) = default;
    ~Quantity() override = default;
    
    using Base::operator==;
    using Base::operator!=;
    bool operator==(double x) const;
    bool operator!=(double x) const;
    
    operator double() const;
};

template<typename T>
struct QuantityContainerTrait<
        T,
        typename std::enable_if<std::is_arithmetic_v<T>>::type
    >
{
    using Type = Quantity;
};

// Quantity holds a scalar value: its common type is the other one
template<typename T>
struct CommonQuantityTypeTrait<
    T, Quantity, 
    // Disable <Quantity, Quantity> specialization to avoid ambiguity
    typename std::enable_if<!std::is_same_v<T, Quantity>>::type>
{
    using Type = T;
};

// Same as above
template<typename T>
struct CommonQuantityTypeTrait<
    Quantity, T,
    typename std::enable_if<!std::is_same_v<T, Quantity>>::type>
{
    using Type = T;
};


/// @brief Compare the magnitude of two compatible quantities
bool operator<(Quantity const & left, Quantity const & right);

/// @brief Compare the magnitude of two compatible quantities
bool operator<=(Quantity const & left, Quantity const & right);

/// @brief Compare the magnitude of two compatible quantities
bool operator>(Quantity const & left, Quantity const & right);

/// @brief Compare the magnitude of two compatible quantities
bool operator>=(Quantity const & left, Quantity const & right);

/// @brief Helper functions for quantity container constructors
namespace details
{

template<std::size_t D>
Dimensions get_dimensions(xt::nested_initializer_list_t<Quantity, D> const & t)
{
    return get_dimensions<D-1>(*t.begin());
}

template<>
Dimensions get_dimensions<1>(xt::nested_initializer_list_t<Quantity, 1> const & t);

template<std::size_t D>
void nested_check_dimensions(
    xt::nested_initializer_list_t<Quantity, D> const & t,
    Dimensions const & dimensions)
{
    return nested_check_dimensions<D-1>(*t.begin(), dimensions);
}

template<>
void nested_check_dimensions<1>(
    xt::nested_initializer_list_t<Quantity, 1> const & t,
    Dimensions const & dimensions);

template<typename Destination>
void nested_copy_magnitude(Destination && destination, Quantity const & source)
{
    *destination = source.magnitude;
    ++destination;
}

template<typename Destination, typename Source>
void nested_copy_magnitude(Destination && destination, Source const & source)
{
    for(auto && item: source)
    {
        nested_copy_magnitude(std::forward<Destination>(destination), item);
    }
}

}

}

namespace std
{

template<>
struct hash<sycomore::Quantity>
{
    std::size_t operator()(sycomore::Quantity const & q) const noexcept;
};

}

// NOTE: forbid Quantity to participate in the operators overloads of namespace xt
namespace xt
{
namespace detail
{
template<typename F, typename T>
struct xfunction_type<F, T, sycomore::Quantity> { };
template<typename F, typename T>
struct xfunction_type<F, sycomore::Quantity, T> { };
template<typename F, typename T>
struct xfunction_type<F, T, sycomore::Quantity const &> { };
template<typename F, typename T>
struct xfunction_type<F, sycomore::Quantity const &, T> { };
}
}

#endif // _dfbc0517_611a_4989_a51c_fa60b94c587f
