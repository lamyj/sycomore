#ifndef _fbf381fe_fd75_427e_88de_a033418c943c
#define _fbf381fe_fd75_427e_88de_a033418c943c

#include <vector>

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
 * @brief Regular EPG model, where the gradient moment is assumed to be
 * identical during each time interval.
 *
 * In this model, the orders of the model are consecutive positive integers
 * starting at 0.
 */
class SYCOMORE_API Regular
{
public:
    Species species;
    
    Regular(
        Species const & species,
        Magnetization const & initial_magnetization={0,0,1}, 
        unsigned int initial_size=100);
    
    Regular(Regular const &) = default;
    Regular(Regular &&) = default;
    Regular & operator=(Regular const &) = default;
    Regular & operator=(Regular &&) = default;
    ~Regular() = default;
    
    /// @brief Return the number of states in the model.
    std::size_t const states_count() const;

    /// @brief Return a given state of the model.
    std::vector<Complex> state(std::size_t order) const;

    /**
     * @brief Return all states in the model, where each state is stored as
     * F̃_k, F̃^*_{-k}, Z̃_k, in order of increasing order.
     */
    std::vector<Complex> states() const;

    /// @brief Return the echo signal, i.e. F̃_0
    Complex const & echo() const;
    
    /// @brief Apply an RF hard pulse.
    void apply_pulse(Quantity angle, Quantity phase=0*units::rad);

    /// @brief Apply a time interval, i.e. relaxation, diffusion, and gradient.
    void apply_time_interval(
        Quantity const & duration, 
        Quantity const & gradient=0*units::T/units::m);

    /// @brief Apply a gradient; in regular EPG, this shifts all orders by 1.
    void shift();

    /// @brief Simulate the relaxation during given duration.
    void relaxation(Quantity const & duration);

    /**
     * @brief Simulate diffusion during given duration with given gradient
     * amplitude.
     */
    void diffusion(Quantity const & duration, Quantity const & gradient);
private:
    std::vector<Complex> _states;
    unsigned int _states_count;
};

}

}

#endif // _fbf381fe_fd75_427e_88de_a033418c943c
