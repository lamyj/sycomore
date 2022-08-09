#ifndef _ce9e6f38_ba33_45d6_8e2a_27552549c830
#define _ce9e6f38_ba33_45d6_8e2a_27552549c830

#include <vector>
#include <xsimd/xsimd.hpp>
#include "sycomore/sycomore.h"

namespace sycomore
{

namespace epg
{

/**
 * @brief Storage classes for single-pool, two pools with exchange, and two
 * pools with magnetization transfer.
 */
namespace pool_storage
{

using Vector = std::vector<Complex, xsimd::aligned_allocator<Complex, 64>>;

/// @brief Single pool: one set of F, F* and Z states.
struct SinglePool
{
    Vector F, F_star, Z;
    
    SinglePool(std::size_t size)
    : F(size), F_star(size), Z(size)
    {
        // Nothing else.
    }
    
    SinglePool(std::size_t size, Complex const & value)
    : F(size, value), F_star(size, value), Z(size, value)
    {
        // Nothing else.
    }
    
    SinglePool(SinglePool const &) = default;
    SinglePool(SinglePool &&) = default;
    SinglePool & operator=(SinglePool const &) = default;
    SinglePool & operator=(SinglePool &&) = default;
    ~SinglePool() = default;
};

/// @brief Two pools with exchange: two separate sets of F, F* and Z states.
struct Exchange
{
    Vector F_a, F_star_a, Z_a, F_b, F_star_b, Z_b;
};

/**
 * @brief Two pools with magnetization transfer: the second pool represents
 * nuclei with very short T2 and thus no transverse magnetization (i.e. no F nor
 * F* states). The longitudinal states are kept separate.
 */
struct MagnetizationTransfer
{
    Vector F, F_star, Z_a, Z_b;
};

}

}

}

#endif // _ce9e6f38_ba33_45d6_8e2a_27552549c830
