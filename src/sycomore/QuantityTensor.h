#ifndef _9b8b809d_7065_4e2a_9c8c_24267f529679
#define _9b8b809d_7065_4e2a_9c8c_24267f529679

#include <xtensor/xtensor.hpp>

#include "sycomore/QuantityBase.h"
// Required for common types
#include "sycomore/QuantityTensorFixed.h"

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
    
    TensorQ(TensorQ<N> const &) = default;
    TensorQ(TensorQ<N> &&) = default;
    TensorQ<N> & operator=(TensorQ<N> const &) = default;
    TensorQ<N> & operator=(TensorQ<N> &&) = default;
    ~TensorQ() override = default;
    
private:
    /// @brief Helper for the constructors
    void _from_array(xt::nested_initializer_list_t<Quantity, N> const & args);
};

template<std::size_t N>
struct QuantityContainerTrait<xt::xtensor<double, N>> { using Type = TensorQ<N>; };

// For TensorFixedQ, the common shape cannot be computed at compile time:
// default to TensorQ with the largest dimension
template<typename S1, typename S2>
struct CommonQuantityTypeTrait<
    TensorFixedQ<S1>, TensorFixedQ<S2>,
    typename std::enable_if<!std::is_same<S1, S2>::value>::type>
{
    using Type = TensorQ<std::max(S1::size(), S2::size())>;
};

// Similar as above
template<typename S, std::size_t N>
struct CommonQuantityTypeTrait<TensorFixedQ<S>, TensorQ<N>>
{
    using Type = TensorQ<std::max(S::size(), N)>;
};

// Similar as above
template<std::size_t N, typename S>
struct CommonQuantityTypeTrait<TensorQ<N>, TensorFixedQ<S>>
{
    using Type = TensorQ<std::max(S::size(), N)>;
};

// Similar as above
template<std::size_t N1, std::size_t N2>
struct CommonQuantityTypeTrait<
    TensorQ<N1>, TensorQ<N2>,
    typename std::enable_if<N1 != N2>::type>
{
    using Type = TensorQ<std::max(N1, N2)>;
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
