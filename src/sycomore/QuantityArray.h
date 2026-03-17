#ifndef _ef70a062_eaab_4ac7_abe8_a22cc86789cd
#define _ef70a062_eaab_4ac7_abe8_a22cc86789cd

#if __has_include(<xtensor/xtensor.hpp>)
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#else
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>
#endif

#include "sycomore/QuantityBase.h"

namespace sycomore
{

/// @brief Quantity container based on xt::xarray
class ArrayQ: public QuantityBase<ArrayQ, xt::xarray<double>>
{
public:
    using Self = ArrayQ;
    using Container = xt::xarray<double>;
    using Base = QuantityBase<Self, Container>;
    using shape_type = Container::shape_type;
    
    using Base::Base;
    
    /// @brief Create an unitialized quantity from a shape
    ArrayQ(shape_type const & shape)
    : Base(Container(shape))
    {
        // Nothing else.
    }
    
    /// @brief Create a scalar quantity
    ArrayQ(xt::nested_initializer_list_t<double, 1> t);
    
    /// @brief Create a quantity container from a homogeneous container of quantity
    ArrayQ(xt::nested_initializer_list_t<Quantity, 1> args);
    
    /// @brief Create a scalar quantity
    ArrayQ(xt::nested_initializer_list_t<double, 2> t);
    
    /// @brief Create a quantity container from a homogeneous container of quantity
    ArrayQ(xt::nested_initializer_list_t<Quantity, 2> args);
    
    /// @brief Create a scalar quantity
    ArrayQ(xt::nested_initializer_list_t<double, 3> t);
    
    /// @brief Create a quantity container from a homogeneous container of quantity
    ArrayQ(xt::nested_initializer_list_t<Quantity, 3> args);
    
    /// @brief Create a scalar quantity
    ArrayQ(xt::nested_initializer_list_t<double, 4> t);
    
    /// @brief Create a quantity container from a homogeneous container of quantity
    ArrayQ(xt::nested_initializer_list_t<Quantity, 4> args);
    
    /// @brief Create a scalar quantity
    ArrayQ(xt::nested_initializer_list_t<double, 5> t);
    
    /// @brief Create a quantity container from a homogeneous container of quantity
    ArrayQ(xt::nested_initializer_list_t<Quantity, 5> args);
    
    ArrayQ(ArrayQ const &) = default;
    ArrayQ(ArrayQ &&) = default;
    ArrayQ & operator=(ArrayQ const &) = default;
    ArrayQ & operator=(ArrayQ &&) = default;
    ~ArrayQ() override = default;
    
private:
    /// @brief Helper for the constructors
    template<std::size_t D>
    void _from_array(xt::nested_initializer_list_t<Quantity, D> const & args);
};

template<typename T>
struct QuantityContainerTrait<
        xt::xarray<T>,
        typename std::enable_if<std::is_arithmetic_v<T>>::type
    >
{
    using Type = ArrayQ;
};

}

namespace std
{

template<>
struct hash<sycomore::ArrayQ>
{
    std::size_t operator()(sycomore::ArrayQ const & q) const noexcept;
};

}

// NOTE: forbid ArrayQ to participate in the operators overloads of namespace xt
namespace xt
{
namespace detail
{
template<typename F, typename T>
struct xfunction_type<F, T, sycomore::ArrayQ> { };
template<typename F, typename T>
struct xfunction_type<F, sycomore::ArrayQ, T> { };
template<typename F, typename T>
struct xfunction_type<F, T, sycomore::ArrayQ const &> { };
template<typename F, typename T>
struct xfunction_type<F, sycomore::ArrayQ const &, T> { };
}
}

#include "QuantityArray.txx"

#endif // _ef70a062_eaab_4ac7_abe8_a22cc86789cd
