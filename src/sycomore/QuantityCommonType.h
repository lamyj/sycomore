#ifndef _bb5d1b6b_eff5_44f4_801a_8ef080d53b93
#define _bb5d1b6b_eff5_44f4_801a_8ef080d53b93

#include <algorithm>
#include <type_traits>

namespace sycomore
{

// Forward declaration of concrete quantity types
class Quantity;
template<typename S> class TensorFixedQ;
template<std::size_t N> class TensorQ;
class ArrayQ;

// Common types to T1 and T2, default to ArrayQ which may hold any dimension
template<typename T1, typename T2, typename Enable = void>
struct CommonQuantityTypeTrait
{
    using Type = ArrayQ;
};

// Helper to pass cv-qualified types down to CommonQuantityTypeTrait
template<typename T1, typename T2>
struct CommonQuantityTypeStruct
{
    using Type = typename CommonQuantityTypeTrait<
            std::remove_cv_t<T1>, std::remove_cv_t<T2>
        >::Type;
};

// User-friendly alias
template<typename T1, typename T2>
using CommonQuantityType = typename CommonQuantityTypeStruct<T1, T2>::Type;

/****************************** Specializations ******************************/

// If the two types are the same, they are their common type
template<typename T>
struct CommonQuantityTypeTrait<T, T>
{
    using Type = T;
};

// Quantity holds a scalar value: its common type is the other one
template<typename T>
struct CommonQuantityTypeTrait<
    T, Quantity, 
    // Disable <Quantity, Quantity> specialization to avoid ambiguity
    typename std::enable_if<!std::is_same<T, Quantity>::value>::type>
{
    using Type = T;
};

// Same as above
template<typename T>
struct CommonQuantityTypeTrait<
    Quantity, T,
    typename std::enable_if<!std::is_same<T, Quantity>::value>::type>
{
    using Type = T;
};

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

#endif // _bb5d1b6b_eff5_44f4_801a_8ef080d53b93
