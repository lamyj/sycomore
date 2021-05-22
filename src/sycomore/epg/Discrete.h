#ifndef _d9169a5f_d53b_4440_bfc7_2b3f978b665d
#define _d9169a5f_d53b_4440_bfc7_2b3f978b665d

#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/epg/robin_hood.h"
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
    Real threshold;
    Quantity delta_omega=0*units::Hz;
    Quantity velocity=0*units::m/units::s;
    
    Discrete(
        Species const & species, 
        Magnetization const & initial_magnetization={0,0,1}, 
        Quantity bin_width=1*units::rad/units::m,
        Real threshold=0);
    
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
     * F(k), F*(-k), Z(k), in order of increasing order.
     */
    std::vector<Complex> states() const;

    /// @brief Return the echo signal, i.e. F_0
    Complex const & echo() const;
    
    /// @brief Apply an RF hard pulse.
    void apply_pulse(Quantity angle, Quantity phase=0*units::rad);

    /** 
     * @brief Apply a time interval, i.e. relaxation, diffusion, gradient, and
     * off-resonance effects.
     */
    void apply_time_interval(
        Quantity const & duration, 
        Quantity const & gradient=0*units::T/units::m, Real threshold=0);
    
    /** 
     * @brief Apply a time interval, i.e. relaxation, diffusion, gradient, and
     * off-resonance effects.
     */
    void apply_time_interval(TimeInterval const & interval);

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
    
    /**
     * @brief Simulate field- and species-related off-resonance effects during 
     * given duration with given frequency offset.
     */
    void off_resonance(Quantity const & duration);
    
    /// @brief Return the bin width.
    Quantity const & bin_width() const;
    
    void bulk_motion(Quantity const & duration, Quantity const & gradient);
    
private:
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> _F, _F_star, _Z;
    Real _M_z_eq;
    
    Quantity _bin_width;
    std::vector<long long, xsimd::aligned_allocator<long long, 64>> _orders;
    
    // Data kept to avoid expansive re-allocation of memory.
    class Cache
    {
    public:
        // Shift-related data.
        // Mapping between a normalized (i.e. folded) order and its location in
        // the states vectors.
        robin_hood::unordered_flat_map<long long, std::size_t> locations;
        decltype(Discrete::_orders) orders;
        decltype(Discrete::_F) F;
        decltype(Discrete::_F_star) F_star;
        decltype(Discrete::_Z) Z;
        
        // Diffusion-related data.
        std::vector<Real, xsimd::aligned_allocator<Real, 64>> k;
        
        void update_shift(std::size_t size);
        void update_diffusion(
            std::size_t size, decltype(Discrete::_orders) const & orders,
            Real bin_width);
        std::size_t get_location(long long order);
    };
    
    Cache _cache;
};
    
}

}

#endif // _d9169a5f_d53b_4440_bfc7_2b3f978b665d
