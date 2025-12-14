#ifndef _55a43ad4_6f2d_4347_8c43_1a9558173c62
#define _55a43ad4_6f2d_4347_8c43_1a9558173c62

#include <xtensor/xview.hpp>
#include <xtensor/xstrided_view.hpp>

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityInterface.h"
#include "sycomore/Quantity.h"

namespace sycomore
{

class Quantity;

/// @brief Read-only view to a quantity magnitude
template<typename V>
class QuantityConstView: public QuantityInterface<QuantityConstView<V>>
{
public:
    using Container = double;
    
    V magnitude;
    Dimensions const & dimensions;
    
    QuantityConstView(V magnitude, Dimensions const & dimensions);
    QuantityConstView(QuantityConstView<V> const & other) = default;
    QuantityConstView(QuantityConstView<V> && other) = default;
    QuantityConstView<V> & operator=(QuantityConstView<V> && other) = delete;
    ~QuantityConstView() override = default;
    
    /// @brief Return the number of elements in a quantity array
    auto size() const;
    
    /// @brief Return the shape of a quantity array
    auto shape() const;
    
    /// @brief Return the d^th dimension of the shape of a quantity array
    auto shape(std::size_t d) const;
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename ... Args>
    QuantityConstReference operator()(Args && ... args) const;
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array, after dimension and bounds checking
     */
    template<typename ... Args>
    QuantityConstReference at(Args && ... args) const;
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename ... Args>
    QuantityConstReference unchecked(Args && ... args) const;
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename Index>
    QuantityConstReference operator[](Index && index) const;
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename T>
    QuantityConstReference operator[](std::initializer_list<T> && index) const;
    
    operator Quantity() const;
    
    template<typename Q>
    operator Q() const;
};

template<typename Q, typename ... S>
auto view(Q const & q, S && ... slices);

template<typename Q>
auto strided_view(Q const & q, xt::xstrided_slice_vector slices);

}

#include "QuantityConstView.txx"

#endif // _55a43ad4_6f2d_4347_8c43_1a9558173c62
