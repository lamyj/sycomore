#ifndef _c4bbf63f_2cdd_4724_9b18_a8dfa6d8fafa
#define _c4bbf63f_2cdd_4724_9b18_a8dfa6d8fafa

#include <type_traits>

namespace sycomore
{

/// @brief Check whether T is a quantity or quantity container type
template<typename T, typename enabled=void>
struct is_quantity: public std::false_type { };

/**
 * @brief Wrapper around std::enable_if for quantity or quantity container types
 */
template<typename T>
using enable_if_quantity = std::enable_if_t<is_quantity<std::decay_t<T>>::value, bool>;

/**
 * @brief Wrapper for std::enable_if for non-quantity and non-quantity container
 * types*/
template<typename T>
using enable_if_not_quantity = std::enable_if_t<!is_quantity<std::decay_t<T>>::value, bool>;

/**
 * @brief Equivalent type which owns its data
 * @sa QuantityReference
 * @sa QuantityConstReference
 */
template<typename T>
struct OwningTypeTrait { using Type = T; };

/// @brief Helper for OwningTypeTrait
template<typename T>
using OwningType = typename OwningTypeTrait<T>::Type;

/// Quantity type associated with container type
template<typename ContainerType>
struct QuantityContainerTrait {};

/// @brief Helper for QuantityContainerTrait
template<typename T>
using QuantityContainer = typename QuantityContainerTrait<T>::Type;

// Forward declaration of concrete quantity type
class ArrayQ;

/// @brief Common types to T1 and T2, default to ArrayQ which may hold any dimension
template<typename T1, typename T2, typename Enable = void>
struct CommonQuantityTypeTrait
{
    using Type = ArrayQ;
};

/// @brief Helper to CommonQuantityTypeTrait for cv-qualified types
template<typename T1, typename T2>
struct CommonQuantityTypeStruct
{
    using Type = typename CommonQuantityTypeTrait<
            OwningType<std::remove_cv_t<T1>>, OwningType<std::remove_cv_t<T2>>
        >::Type;
};

/// @brief Helper for CommonQuantityTypeStruct
template<typename T1, typename T2>
using CommonQuantityType = typename CommonQuantityTypeStruct<T1, T2>::Type;

// If the two types are the same, they are their common type
template<typename T>
struct CommonQuantityTypeTrait<T, T>
{
    using Type = T;
};

}

#endif // _c4bbf63f_2cdd_4724_9b18_a8dfa6d8fafa
