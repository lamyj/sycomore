#ifndef _fbf381fe_fd75_427e_88de_a033418c943c
#define _fbf381fe_fd75_427e_88de_a033418c943c

#include <vector>

#include <xsimd/xsimd.hpp>

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
 * @brief Regular EPG model, where the gradient moment is assumed to be
 * a multiple of a user-specified unitary moment during each time interval.
 *
 * In this model, the orders of the model are consecutive positive integers
 * starting at 0.
 */
class SYCOMORE_API Regular
{
public:
    Species species;
    Real threshold=0;
    Quantity delta_omega=0*units::Hz;
    Quantity velocity=0*units::m/units::s;
    
    Regular(
        Species const & species, 
        Magnetization const & initial_magnetization={0,0,1}, 
        unsigned int initial_size=100, 
        Quantity const & unit_gradient_area=0*units::mT/units::m*units::ms,
        double gradient_tolerance=1e-5);
    
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
     * F_k, F*_{-k}, Z_k, in order of increasing order.
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
        Quantity const & gradient=0*units::T/units::m);
    
    /** 
     * @brief Apply a time interval, i.e. relaxation, diffusion, gradient, and
     * off-resonance effects.
     */
    void apply_time_interval(TimeInterval const & interval);

    /// @brief Apply a unit gradient; in regular EPG, this shifts all orders by 1.
    void shift();
    
    /* 
     * @brief Apply an arbitrary gradient; in regular EPG, this shifts all 
     * orders by an integer number corresponding to a multiple of the unit 
     * gradient.
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
    
    void bulk_motion(Quantity const & duration, Quantity const & gradient);
    
    Quantity const & unit_gradient_area() const;
    double gradient_tolerance() const;
    
private:
    std::vector<Complex, xsimd::aligned_allocator<Complex, 64>> _F, _F_star, _Z;
    Real _M_z_eq;
    unsigned int _states_count;
    
    /// @brief Area of the unit gradient, in T/m*s.
    Quantity _unit_gradient_area;
    
    /** 
     * @brief Tolerance used when checking that a prescribed gradient moment is
     * close enough to a multiple of the unit gradient.
     */
    double _gradient_tolerance;
    
    /// @brief Shift all orders by given number of steps (may be negative).
    void _shift(int n);
    
    // Data kept to avoid expansive re-allocation of memory.
    class Cache
    {
    public:
        // Diffusion-related data.
        std::vector<Real, xsimd::aligned_allocator<Real, 64>> k;
        
        void update_diffusion(std::size_t size, Real unit_gradient_area);
    };
    
    Cache _cache;
};

}

}

#endif // _fbf381fe_fd75_427e_88de_a033418c943c
