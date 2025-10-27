#ifndef _9b8b809d_7065_4e2a_9c8c_24267f529679
#define _9b8b809d_7065_4e2a_9c8c_24267f529679

#include <xtensor/xtensor.hpp>

#include "sycomore/QuantityBase.h"

namespace sycomore
{

/// @brief Quantity container based on xt::xtensor
template<std::size_t N>
class TensorQ: public QuantityBase<TensorQ<N>, xt::xtensor<double, N>>
{
public:
    using Self = TensorQ<N>;
    using Container = xt::xtensor<double, N>;
    using Base = QuantityBase<Self, Container>;
    
    using Base::Base;
    
    /// @brief Create a scalar quantity
    TensorQ(xt::nested_initializer_list_t<double, N> t);
    
    /// @brief Create a quantity container from a homogeneous container of quantity
    TensorQ(xt::nested_initializer_list_t<Quantity, N> args);
    
private:
    /// @brief Helper for the constructors
    void _from_array(xt::nested_initializer_list_t<Quantity, N> const & args);
};

}

namespace std
{

template<std::size_t N>
struct hash<sycomore::TensorQ<N>>
{
    std::size_t operator()(sycomore::TensorQ<N> const & q) const noexcept;
};

}

#include "QuantityTensor.txx"

#endif // _9b8b809d_7065_4e2a_9c8c_24267f529679
