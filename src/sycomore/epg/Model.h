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

/// @brief EPG simulators
namespace epg
{

/// @brief Model a single- or two-pool system expressed in EPG formalism
class Model
{
public:
    /// @brief Populations of the states
    using Population = std::vector<
        Complex, xsimd::aligned_allocator<Complex, 64>>;
    
    /// @brief Kind of the model, determines the number of pools
    enum Kind { SinglePool, Exchange, MagnetizationTransfer };
    
    /// @brief Kind of the model, determines the number of pools
    Kind const kind;
    
    /// @brief Number of pools
    std::size_t const pools;
    
    /// @brief Species
    std::vector<Species> species;
    
    /// @brief Equilibrium magnetization on the z axis
    std::vector<Real> M0;
    
    /// @brief Exchange rates between the pools, in Hz
    std::vector<Quantity> k;
    
    /// @brief Frequency offset of pool b w.r.t. to pool al.
    Quantity delta_b;
    
    /// @brief EPG F states for each pool
    std::vector<Population> F;
    /// @brief EPG F* states for each pool
    std::vector<Population> F_star;
    /// @brief EPG Z states for each pool
    std::vector<Population> Z;
    
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
    
    /// @brief Default copy constructor
    Model(Model const &) = default;
    /// @brief Default move constructor
    Model(Model &&) = default;
    /// @brief Default copy assignment
    Model & operator=(Model const & other);
    /// @brief Default move assignment
    Model & operator=(Model && other);
    
private:
    void _initialize(
        std::vector<Vector3R> const & M0, std::size_t initial_size);
};
    
}

}

#endif // _bbf58904_14a3_486b_8a48_7c982bcb7955
