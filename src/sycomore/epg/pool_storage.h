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
class Base
{
public:
    Vector F, F_star, Z;
    
    Base(std::size_t size);
    Base(std::size_t size, Complex const & value);
    
    Base(Base const &) = default;
    Base(Base &&) = default;
    Base & operator=(Base const &) = default;
    Base & operator=(Base &&) = default;
    virtual ~Base() = default;
};

using SinglePool = Base;

/// @brief Two pools with exchange: two separate sets of F, F* and Z states.
class Exchange: public Base
{
public:
    Vector & F_a, & F_star_a, &Z_a;
    Vector F_b, F_star_b, Z_b;
    
    Exchange(std::size_t size);
    Exchange(std::size_t size, Complex const & value);
    
    Exchange(Exchange const & other);
    Exchange(Exchange && other);
    Exchange & operator=(Exchange const & other);
    Exchange & operator=(Exchange && other);
    
    virtual ~Exchange() = default;
};

/**
 * @brief Two pools with magnetization transfer: the second pool represents
 * nuclei with very short T2 and thus no transverse magnetization (i.e. no F nor
 * F* states). The longitudinal states are kept separate.
 */
struct MagnetizationTransfer: public Base
{
    Vector & Z_a;
    Vector Z_b;
    
    MagnetizationTransfer(std::size_t size);
    MagnetizationTransfer(std::size_t size, Complex const & value);
    
    MagnetizationTransfer(MagnetizationTransfer const & other);
    MagnetizationTransfer(MagnetizationTransfer && other);
    MagnetizationTransfer & operator=(MagnetizationTransfer const & other);
    MagnetizationTransfer & operator=(MagnetizationTransfer && other);
    
    virtual ~MagnetizationTransfer() = default;
};

}

}

}

#endif // _ce9e6f38_ba33_45d6_8e2a_27552549c830
