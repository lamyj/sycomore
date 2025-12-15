#ifndef _9b8b809d_7065_4e2a_9c8c_24267f529679
#define _9b8b809d_7065_4e2a_9c8c_24267f529679

#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>

#include "sycomore/QuantityArray.h"
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
    using shape_type = typename Container::shape_type;
    
    using Base::Base;
    
    /// @brief Create an unitialized quantity from a shape
    TensorQ(shape_type const & shape)
    : Base(Container(shape))
    {
        // Nothing else.
    }
    
    /// @brief Create a quantity from a shape-compatible Quantity array
    TensorQ(ArrayQ const & other);
    
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

template<typename T, std::size_t N>
struct QuantityContainerTrait<
        xt::xtensor<T, N>,
        typename std::enable_if<std::is_arithmetic_v<T>>::type>
{
    using Type = TensorQ<N>;
};

// For TensorFixedQ, the common shape cannot be computed at compile time:
// default to TensorQ with the largest dimension
template<typename S1, typename S2>
struct CommonQuantityTypeTrait<
    TensorFixedQ<S1>, TensorFixedQ<S2>,
    typename std::enable_if<!std::is_same_v<S1, S2>>::type>
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

// NOTE: forbid TensorQ to participate in the operators overloads of namespace xt
namespace xt
{
namespace detail
{
template<typename F, typename T, std::size_t N>
struct xfunction_type<F, T, sycomore::TensorQ<N>> { };
template<typename F, std::size_t N, typename T>
struct xfunction_type<F, sycomore::TensorQ<N>, T> { };
template<typename F, typename T, std::size_t N>
struct xfunction_type<F, T, sycomore::TensorQ<N> const &> { };
template<typename F, std::size_t N, typename T>
struct xfunction_type<F, sycomore::TensorQ<N> const &, T> { };
}
}


#include "QuantityTensor.txx"

#endif // _9b8b809d_7065_4e2a_9c8c_24267f529679
