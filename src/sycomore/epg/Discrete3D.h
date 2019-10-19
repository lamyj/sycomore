#ifndef _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057
#define _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057

#include <vector>

#include "sycomore/Array.h"
#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

/**
 * @brief Discrete EPG in three dimensions.
 *
 * @note To keep the "folded" aspect of 1D EPG (i.e. the states being stored as
 * (F̃_{+}(-k), F̃^{*}_{-}(-k), Z̃(k)) with k>=), the 3D EPG must store the states
 * of a half-space bounded by an arbitrary plane including the origin. This
 * arbitrary plane is chosen as the x=0 plane.
 */
class SYCOMORE_API Discrete3D
{
public:
    using Order = Array<Quantity>;
    using State = std::vector<Complex>;

    Species species;

    Discrete3D(
        Species const & species,
        Magnetization const & initial_magnetization={0,0,1},
        Quantity bin_width=1*units::rad/units::m);

    Discrete3D(Discrete3D const &) = default;
    Discrete3D(Discrete3D &&) = default;
    Discrete3D & operator=(Discrete3D const &) = default;
    Discrete3D & operator=(Discrete3D &&) = default;
    ~Discrete3D() = default;

    /// @brief Return the number of states of the model.
    std::size_t size() const;

    /// @brief Return the orders of the model.
    std::vector<Quantity> orders() const;

    /// @brief Return a given state of the model.
    State state(Order const & order) const;

    /**
     * @brief Return all states in the model, where each state is stored as
     * F̃(k), Z̃(k).
     */
    std::vector<Complex> const & states() const;

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
    
    /// @brief Return the bin width.
    Quantity const & bin_width() const;

private:
    using Bin = Array<int64_t>;
    std::vector<Bin::value_type> _orders;
    std::vector<Complex> _states;
    decltype(_states)::iterator _zero_it;
    Quantity _bin_width;
};

}

}

#endif // _fcca9c67_7c2f_4a9d_abbb_718dc5fd0057
