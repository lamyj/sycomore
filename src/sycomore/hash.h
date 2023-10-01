#ifndef _a26d369d_eae0_467a_98b4_dde5e537b8ec
#define _a26d369d_eae0_467a_98b4_dde5e537b8ec

#include <cstddef>
#include <functional>

namespace sycomore
{

/// @brief Combine two hashes, implementation from boost::hash_combine.
void combine_hashes(std::size_t & seed, std::size_t value);

}

#endif // _a26d369d_eae0_467a_98b4_dde5e537b8ec
