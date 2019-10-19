#ifndef _d9169a5f_d53b_4440_bfc7_2b3f978b665d
#define _d9169a5f_d53b_4440_bfc7_2b3f978b665d

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
 * @brief Discrete EPG model, where the gradient moments may vary across time
 * intervals.
 *
 * In this model, the orders of the model are stored in bins of user-specified
 * width (hence the term "discrete"), expressed in rad/m.
 */
class SYCOMORE_API Discrete
{
public:
    using Order = Quantity;
    using State = std::vector<Complex>;

    Species species;
    
    Discrete(
        Species const & species, 
        Magnetization const & initial_magnetization={0,0,1}, 
        Quantity bin_width=1*units::rad/units::m);
    
    Discrete(Discrete const &) = default;
    Discrete(Discrete &&) = default;
    Discrete & operator=(Discrete const &) = default;
    Discrete & operator=(Discrete &&) = default;
    ~Discrete() = default;

    /// @brief Return the number of states of the model.
    std::size_t size() const;
    
    /// @brief Return the orders of the model.
    std::vector<Quantity> orders() const;

    /// @brief Return a given state of the model.
    std::vector<Complex> state(std::size_t bin) const;

    /// @brief Return a given state of the model.
    std::vector<Complex> state(Quantity const & order) const;

    /**
     * @brief Return all states in the model, where each state is stored as
     * F̃(k), F̃^*(-k), Z̃(k), in order of increasing order.
     */
    std::vector<Complex> const & states() const;

    /// @brief Return the echo signal, i.e. F̃_0
    Complex const & echo() const;
    
    /// @brief Apply an RF hard pulse.
    void apply_pulse(Quantity angle, Quantity phase=0*units::rad);

    /// @brief Apply a time interval, i.e. relaxation, diffusion, and gradient.
    void apply_time_interval(
        Quantity const & duration, 
        Quantity const & gradient=0*units::T/units::m, Real threshold=0);

    /**
     * @brief Apply a gradient; in discrete EPG, this shifts all orders by
     * specified value.
     */
    void shift(Quantity const & duration, Quantity const & gradient);

    /// @brief Simulate the relaxation during given duration.
    void relaxation(Quantity const & duration);

    /**
     * @brief Simulate diffusion during given duration with given gradient
     * amplitude.
     */
    void diffusion(Quantity const & duration, Quantity const & gradient);
    
    /// @brief Return the bin width.
    Quantity const & bin_width() const;

private:
    std::vector<Complex> _states;
    Quantity _bin_width;
    std::vector<long long> _orders;
};
    
}

}

#endif // _d9169a5f_d53b_4440_bfc7_2b3f978b665d
