#ifndef _bbf58904_14a3_486b_8a48_7c982bcb7955
#define _bbf58904_14a3_486b_8a48_7c982bcb7955

#include <vector>
#include <xsimd/xsimd.hpp>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"

namespace sycomore
{

namespace epg
{

/// @brief Model a single- or two-pool system expressed in EPG formalism
class Model
{
public:
    using Population = std::vector<
        Complex, xsimd::aligned_allocator<Complex, 64>>;
    enum Kind { SinglePool, Exchange, MagnetizationTransfer };
    
    Kind const kind;
    std::size_t const pools;
    
    /// @brief
    std::vector<Species> species;
    
    /// @brief Equilibrium magnetization on the z axis
    std::vector<Real> M0;
    
    /// @brief Exchange rates between the pools, in Hz
    std::vector<Quantity> k;
    
    /// @brief Frequency offset of pool b w.r.t. to pool al.
    Quantity delta_b;
    
    /// @brief EPG states for each pool
    std::vector<Population> F, F_star, Z;
    
    /// @brief Create a single-pool model.
    Model(
        Species const & species, Vector3R const & M0,
        std::size_t initial_size);
    
    /// @brief Create an exchange model.
    Model(
        Species const & species_a, Species const & species_b,
        Vector3R const & M0_a, Vector3R const & M0_b,
        Quantity const & k_a, Quantity const & delta_b,
        std::size_t initial_size);
    
    /// @brief Create a magnetization transfer model.
    Model(
        Species const & species_a, Quantity const & R1_b_or_T1_b,
        Vector3R const & M0_a, Vector3R const & M0_b,
        Quantity const & k_a,
        std::size_t initial_size);
    
    Model(Model const &) = default;
    Model(Model &&) = default;
    Model & operator=(Model const & other);
    Model & operator=(Model && other);
    
private:
    void _initialize(
        std::vector<Vector3R> const & M0, std::size_t initial_size);
};
    
}

}

#endif // _bbf58904_14a3_486b_8a48_7c982bcb7955
