#ifndef _d9169a5f_d53b_4440_bfc7_2b3f978b665d
#define _d9169a5f_d53b_4440_bfc7_2b3f978b665d

#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/epg/Base.h"
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
class SYCOMORE_API Discrete: public Base
{
public:
    using Order = Quantity;
    using State = std::vector<Complex>;
    
    Quantity velocity=0*units::m/units::s;
    
    Discrete(
        Species const & species, 
        Magnetization const & initial_magnetization={0,0,1}, 
        Quantity bin_width=1*units::rad/units::m);
    
    Discrete(Discrete const &) = default;
    Discrete(Discrete &&) = default;
    Discrete & operator=(Discrete const &) = default;
    Discrete & operator=(Discrete &&) = default;
    virtual ~Discrete() = default;

    /// @brief Return the number of states of the model.
    virtual std::size_t size() const;
    
    /// @brief Return the orders of the model.
    std::vector<Order> orders() const;
    
    using Base::state;
    
    /// @brief Return a given state of the model.
    std::vector<Complex> state(Order const & order) const;

    /** 
     * @brief Apply a time interval, i.e. relaxation, diffusion, gradient, and
     * off-resonance effects.
     */
    void apply_time_interval(
        Quantity const & duration, 
        Quantity const & gradient=0*units::T/units::m);
    
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

    /**
     * @brief Simulate diffusion during given duration with given gradient
     * amplitude.
     */
    void diffusion(Quantity const & duration, Quantity const & gradient);
    
    /// @brief Return the bin width.
    Quantity const & bin_width() const;
    
    void bulk_motion(Quantity const & duration, Quantity const & gradient);
    
private:
    using Orders = std::vector<long long, xsimd::aligned_allocator<long long, 64>>;
    Quantity _bin_width;
    Orders _orders;
    
    // Data kept to avoid expansive re-allocation of memory.
    class Cache
    {
    public:
        // Shift-related data.
        // Mapping between a normalized (i.e. folded) order and its location in
        // the states vectors.
        robin_hood::unordered_flat_map<long long, std::size_t> locations;
        Orders orders;
        std::vector<Model::Population> F, F_star, Z;
        
        // Diffusion-related data.
        std::vector<Real, xsimd::aligned_allocator<Real, 64>> k;
        
        Cache(std::size_t pools);
        
        void update_shift(std::size_t size);
        void update_diffusion(
            std::size_t size, Orders const & orders, Real bin_width);
        std::size_t get_location(long long order);
    };
    
    Cache _cache;
};
    
}

}

#endif // _d9169a5f_d53b_4440_bfc7_2b3f978b665d
