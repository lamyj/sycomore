#ifndef _a26d369d_eae0_467a_98b4_dde5e537b8ec
#define _a26d369d_eae0_467a_98b4_dde5e537b8ec

#include <cstddef>
#include <functional>

#include "sycomore/sycomore_api.h"

namespace sycomore
{

/// @brief Combine two hashes, implementation from boost::hash_combine.
SYCOMORE_API void combine_hashes(std::size_t & seed, std::size_t value);

template<typename T>
std::size_t hash_range(T const * begin, T const * end)
{
    std::size_t seed=0;
    std::hash<T> hasher;
    while(begin != end)
    {
        combine_hashes(seed, hasher(*begin));
        ++begin;
    }
    return seed;
}

}

#endif // _a26d369d_eae0_467a_98b4_dde5e537b8ec
