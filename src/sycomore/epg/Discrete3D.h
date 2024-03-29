#ifndef _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057
#define _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057

#include <array>
#include <vector>

#include <xsimd/xsimd.hpp>
#include <xtensor/xarray.hpp>

#include "sycomore/Array.h"
#include "sycomore/Buffer.h"
#include "sycomore/epg/Base.h"
#include "sycomore/epg/robin_hood.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

namespace std
{

/// @brief Hash functor
template<>
struct hash<array<int64_t, 3>>
{
    /// @brief Hash function
    typename std::enable_if<sizeof(std::size_t)==sizeof(int64_t), std::size_t>::type
    operator()(array<int64_t, 3> const & a) const
    {
        std::size_t value = 0;
        value ^= reinterpret_cast<std::size_t const &>(a[0]) + 0x9e3779b9 + (value<<6) + (value>>2);
        value ^= reinterpret_cast<std::size_t const &>(a[1]) + 0x9e3779b9 + (value<<6) + (value>>2);
        value ^= reinterpret_cast<std::size_t const &>(a[2]) + 0x9e3779b9 + (value<<6) + (value>>2);
        return value;
    }
};

}

namespace sycomore
{

namespace epg
{

/**
 * @brief Discrete EPG in which the gradients may be specified in three 
 * dimensions.
 */
class Discrete3D: public Base
{
public:
    /// @brief Order of the model, as gradient area
    using Order = Vector3Q;
    
    /// @brief Create a single-pool model
    Discrete3D(
        Species const & species,
        Vector3R const & initial_magnetization={0,0,1},
        Quantity bin_width=1*units::rad/units::m);
    
    /// @brief Create an exchange model
    Discrete3D(
        Species const & species_a, Species const & species_b,
        Vector3R const & M0_a, Vector3R const & M0_b,
        Quantity const & k_a, Quantity const & delta_b=0*units::Hz,
        Quantity bin_width=1*units::rad/units::m);
    
    /// @brief Create an MT model
    Discrete3D(
        Species const & species_a, Quantity const & R1_b_or_T1_b,
        Vector3R const & M0_a, Vector3R const & M0_b,
        Quantity const & k_a,
        Quantity bin_width=1*units::rad/units::m);

    /// @brief Default copy constructor
    Discrete3D(Discrete3D const &) = default;
    /// @brief Default move constructor
    Discrete3D(Discrete3D &&) = default;
    /// @brief Default copy assignment
    Discrete3D & operator=(Discrete3D const &) = default;
    /// @brief Default move assignment
    Discrete3D & operator=(Discrete3D &&) = default;
    /// @brief Default destructor
    virtual ~Discrete3D() = default;

    /// @brief Return the number of states of the model.
    virtual std::size_t size() const;

    /// @brief Return the orders of the model.
    TensorQ<2> orders() const;
    
    using Base::state;
    
    /// @brief Return a given state of the model.
    ArrayC state(Order const & order) const;

    /// @brief Apply a time interval, i.e. relaxation, diffusion, and gradient.
    void apply_time_interval(
        Quantity const & duration,
        Vector3Q const & gradient={
            0*units::T/units::m,0*units::T/units::m,0*units::T/units::m});
    
    /// @brief Apply a time interval, i.e. relaxation, diffusion, and gradient.
    void apply_time_interval(TimeInterval const & interval);

    /**
     * @brief Apply a gradient; in discrete EPG, this shifts all orders by
     * specified value.
     */
    void shift(Quantity const & duration, Vector3Q const & gradient);

    /**
     * @brief Simulate diffusion during given duration with given gradient
     * amplitude.
     */
    void diffusion(Quantity const & duration, Vector3Q const & gradient);
    
    /// @brief Return the bin width.
    Quantity const & bin_width() const;

private:
    using Bin = std::array<int64_t, 3>;
    using Orders = Buffer<Bin::value_type>;
    Orders _orders;

    Quantity _bin_width;
    
    // Data kept to avoid expansive re-allocation of memory.
    class Cache
    {
    public:
        using RealVector = Buffer<Real>;
        
        // Shift-related data.
        // Mapping between a normalized (i.e. folded) order and its location in
        // the states vectors.
        robin_hood::unordered_flat_map<Bin, std::size_t> locations;
        Orders orders;
        std::vector<Model::Population> F, F_star, Z;
        
        // Diffusion-related data.
        std::vector<RealVector, xsimd::aligned_allocator<RealVector, 64>> k;
        RealVector b_L_D;
        RealVector b_T_plus_D;
        RealVector b_T_minus_D;
        
        Cache(std::size_t pools);
        
        void update_shift(std::size_t size);
        void update_diffusion(
            std::size_t size, Orders const & orders, Real bin_width);
        std::size_t location(Bin const & order);
    };
    
    Cache _cache;
};

}

}

#endif // _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057
