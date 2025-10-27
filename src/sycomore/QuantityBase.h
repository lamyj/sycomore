#ifndef _2bdd9ca4_0f72_4689_927a_df66e17a1a31
#define _2bdd9ca4_0f72_4689_927a_df66e17a1a31

#include <initializer_list>

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityConstReference.h"
#include "sycomore/QuantityInterface.h"
#include "sycomore/QuantityReference.h"

namespace sycomore
{

/// @brief Base class for quantity-like objects which own their data
template<typename TDerived, typename TContainer>
class QuantityBase: public QuantityInterface<TDerived>
{
public:
    using Container = TContainer;
    
    /// @brief Magnitude of the quantity in base units
    Container magnitude;
    
    /// @brief Dimensions of the quantity
    Dimensions dimensions;
    
    /// @brief Create a quantity from a magnitude and dimensions
    QuantityBase(Container const & magnitude={}, Dimensions const & dimensions={});
    
    /// @brief Return the number of elements in a quantity array
    auto size() const
    {
        return this->magnitude.size();
    }
    
    /// @brief Return the shape of a quantity array
    auto shape() const
    {
        return this->magnitude.shape();
    }
    
    /// @brief Return the d^th dimension of the shape of a quantity array
    auto shape(std::size_t d) const
    {
        return this->magnitude.shape(d);
    }
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename ... Args>
    QuantityConstReference operator()(Args && ... args) const;
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename ... Args>
    QuantityReference operator()(Args && ... args);
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array, after dimension and bounds checking
     */
    template<typename ... Args>
    QuantityConstReference at(Args && ... args) const;
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array, after dimension and bounds checking
     */
    template<typename ... Args>
    QuantityReference at(Args && ... args);
    
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
    template<typename ... Args>
    QuantityReference unchecked(Args && ... args);
    
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
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename Index>
    QuantityReference operator[](Index && index);
    
    /**
     * @brief Returns a reference to the element at the specified position of
     * a quantity array
     */
    template<typename T>
    QuantityReference operator[](std::initializer_list<T> && index);
    
    /**
     * @brief Return the scalar value of the quantity converted to the given 
     * unit.
     *
     * Raise an exception if the given unit is not compatible.
     */
    template<typename D, typename C, std::enable_if_t<std::is_scalar<C>::value, bool> = true>
    TContainer convert_to(QuantityBase<D, C> const & destination) const
    {
        this->check_dimensions(destination, "Conversion requires same dimensions");
        return this->magnitude/destination.magnitude;
    }
    
    /**
     * @brief Convert to a scalar.
     *
     * Raise an exception if the quantity is not unitless.
     */
    TContainer const & scalar() const;
};

}

#include "sycomore/QuantityBase.txx"

#endif // _2bdd9ca4_0f72_4689_927a_df66e17a1a31
