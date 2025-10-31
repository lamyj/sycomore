#ifndef _0b6aefc6_cf98_4fff_965c_d22c09e27aac
#define _0b6aefc6_cf98_4fff_965c_d22c09e27aac

#include <xtensor/xfixed.hpp>
#include <xtensor/xio.hpp>

#include "sycomore/Quantity.h"
#include "sycomore/QuantityBase.h"

namespace sycomore
{

    /// @brief Quantity container based on xt::xtensor_fixed
template<typename S>
class TensorFixedQ: public QuantityBase<TensorFixedQ<S>, xt::xtensor_fixed<double, S>>
{
public:
    using Self = TensorFixedQ<S>;
    using Container = xt::xtensor_fixed<double, S>;
    using Base = QuantityBase<Self, Container>;
    using shape_type = typename Container::shape_type;
    
    static constexpr std::size_t const rank = std::tuple_size<S>::value;
    
    using Base::Base;
    
    /// @brief Create a scalar quantity
    TensorFixedQ(xt::nested_initializer_list_t<double, std::tuple_size<S>::value> t);
    
    /// @brief Create a quantity container from a homogeneous container of quantity
    TensorFixedQ(xt::nested_initializer_list_t<Quantity, rank> args);
    
    TensorFixedQ(TensorFixedQ<S> const &) = default;
    TensorFixedQ(TensorFixedQ<S> &&) = default;
    TensorFixedQ<S> & operator=(TensorFixedQ<S> const &) = default;
    TensorFixedQ<S> & operator=(TensorFixedQ<S> &&) = default;
    ~TensorFixedQ() override = default;
private:
    /// @brief Helper for the constructors
    void _from_array(xt::nested_initializer_list_t<Quantity, rank> const & args);
};

template<typename S>
struct QuantityContainerTrait<xt::xtensor_fixed<double, S>>
{
    using Type = TensorFixedQ<S>;
};

using Vector2Q = TensorFixedQ<xt::xshape<2>>;
using Vector3Q = TensorFixedQ<xt::xshape<3>>;
using Vector4Q = TensorFixedQ<xt::xshape<4>>;

using Matrix2x2Q = TensorFixedQ<xt::xshape<2, 2>>;
using Matrix3x3Q = TensorFixedQ<xt::xshape<3, 3>>;

}

namespace std
{

template<typename S>
struct hash<sycomore::TensorFixedQ<S>>
{
    std::size_t operator()(sycomore::TensorFixedQ<S> const & q) const noexcept;
};

}

#include "QuantityTensorFixed.txx"

#endif // _0b6aefc6_cf98_4fff_965c_d22c09e27aac
