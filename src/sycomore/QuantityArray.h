#ifndef _ef70a062_eaab_4ac7_abe8_a22cc86789cd
#define _ef70a062_eaab_4ac7_abe8_a22cc86789cd

#include <xtensor/xarray.hpp>

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
    
    using Base::Base;
    
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
    
private:
    /// @brief Helper for the constructors
    template<std::size_t D>
    void _from_array(xt::nested_initializer_list_t<Quantity, D> const & args);
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

#include "QuantityArray.txx"

#endif // _ef70a062_eaab_4ac7_abe8_a22cc86789cd
