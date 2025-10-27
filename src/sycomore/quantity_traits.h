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
using enable_if_quantity = std::enable_if_t<is_quantity<T>::value, bool>;

/**
 * @brief Wrapper for std::enable_if for non-quantity and non-quantity container
 * types*/
template<typename T>
using enable_if_not_quantity = std::enable_if_t<!is_quantity<T>::value, bool>;

/**
 * @brief Return an equivalent type which owns its data
 * @sa QuantityReference
 * @sa QuantityConstReference
 */
template<typename T>
struct OwningTypeTrait { using Type = T; };

template<typename T>
using OwningType = typename OwningTypeTrait<T>::Type;

}

#endif // _c4bbf63f_2cdd_4724_9b18_a8dfa6d8fafa
