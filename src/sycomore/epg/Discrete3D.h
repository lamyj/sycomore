#ifndef _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057
#define _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057

#include <array>
#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/Array.h"
#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

/**
 * @brief Discrete EPG in which the gradients may be specified in three 
 * dimensions.
 */
class SYCOMORE_API Discrete3D
{
public:
    using Order = Array<Quantity>;
    using State = std::vector<Complex>;

    Species species;
    Real threshold;
    Quantity delta_omega=0*units::Hz;

    Discrete3D(
        Species const & species,
        Magnetization const & initial_magnetization={0,0,1},
        Quantity bin_width=1*units::rad/units::m,
        Real threshold=0);

    Discrete3D(Discrete3D const &) = default;
    Discrete3D(Discrete3D &&) = default;
    Discrete3D & operator=(Discrete3D const &) = default;
    Discrete3D & operator=(Discrete3D &&) = default;
    ~Discrete3D() = default;

    /// @brief Return the number of states of the model.
    std::size_t size() const;

    /// @brief Return the orders of the model.
    std::vector<Quantity> orders() const;
    
    State state(std::size_t order) const;

    /// @brief Return a given state of the model.
    State state(Order const & order) const;

    /**
     * @brief Return all states in the model, where each state is stored as
     * F̃(k), Z̃(k).
     */
    std::vector<Complex> states() const;

    /// @brief Return the echo signal, i.e. F̃_0
    Complex const & echo() const;

    /// @brief Apply an RF hard pulse.
    void apply_pulse(Quantity angle, Quantity phase=0*units::rad);

    /// @brief Apply a time interval, i.e. relaxation, diffusion, and gradient.
    void apply_time_interval(
        Quantity const & duration,
        Array<Quantity> const & gradient={
            0*units::T/units::m,0*units::T/units::m,0*units::T/units::m,},
        Real threshold=0);
    
    /// @brief Apply a time interval, i.e. relaxation, diffusion, and gradient.
    void apply_time_interval(TimeInterval const & interval);

    /**
     * @brief Apply a gradient; in discrete EPG, this shifts all orders by
     * specified value.
     */
    void shift(Quantity const & duration, Array<Quantity> const & gradient);

    /// @brief Simulate the relaxation during given duration.
    void relaxation(Quantity const & duration);

    /**
     * @brief Simulate diffusion during given duration with given gradient
     * amplitude.
     */
    void diffusion(Quantity const & duration, Array<Quantity> const & gradient);
    
    /**
     * @brief Simulate field- and species-related off-resonance effects during 
     * given duration with given frequency offset.
     */
    void off_resonance(Quantity const & duration);
    
    /// @brief Return the bin width.
    Quantity const & bin_width() const;

private:
    using Bin = std::array<int64_t, 3>;
    std::vector<Bin::value_type, xsimd::aligned_allocator<Bin::value_type, 64>> _orders;
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> _F;
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> _F_star;
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> _Z;

    Quantity _bin_width;
};

}

}

#endif // _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057
