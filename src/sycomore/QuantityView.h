#ifndef _c86da767_3b1b_44e0_ad8e_6c294a8e4df9
#define _c86da767_3b1b_44e0_ad8e_6c294a8e4df9

#if __has_include(<xtensor/xtensor.hpp>)
#include <xtensor/xview.hpp>
#else
#include <xtensor/views/xview.hpp>
#endif

#include "sycomore/QuantityBase.h"

namespace sycomore
{

/// @brief View to a part of a quantity container
template<typename QuantityContainer, typename View>
class QuantityView: public QuantityBase<QuantityView<QuantityContainer, View>, View>
{
public:
    using Self = QuantityView<QuantityContainer, View>;
    using Container = View;
    using Base = QuantityBase<Self, Container>;
    
    using Base::Base;
    
    /// @brief Convert to a concrete container
    template<typename T, enable_if_quantity<T> = true>
    operator T() const
    {
        return T{this->magnitude, this->dimensions};
    }
};

/// @brief Return a view to a part of a quantity container
template<typename QuantityContainer, typename ... Slices>
auto
view(QuantityContainer && container, Slices && ... slices)
{
    auto magnitude = xt::view(container.magnitude, slices ...);
    return QuantityView<
            QuantityContainer, decltype(magnitude)
        >(magnitude, container.dimensions);
}

}

#endif // _c86da767_3b1b_44e0_ad8e_6c294a8e4df9
