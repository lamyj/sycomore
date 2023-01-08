#ifndef _fbf381fe_fd75_427e_88de_a033418c943c
#define _fbf381fe_fd75_427e_88de_a033418c943c

#include <vector>

#include <xsimd/xsimd.hpp>

#include "sycomore/Array.h"
#include "sycomore/epg/Base.h"
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
class SYCOMORE_API Regular: public Base
{
public:
    using Order = Quantity;
    
    Quantity velocity=0*units::m/units::s;
    
    Regular(
        Species const & species, 
        Vector3<Real> const & initial_magnetization={0,0,1}, 
        unsigned int initial_size=100, 
        Quantity const & unit_gradient_area=0*units::mT/units::m*units::ms,
        double gradient_tolerance=1e-5);
    
    Regular(
        Species const & species_a, Species const & species_b,
        Vector3<Real> const & M0_a, Vector3<Real> const & M0_b,
        Quantity const & k_a, Quantity const & delta_b=0*units::Hz,
        unsigned int initial_size=100, 
        Quantity const & unit_gradient_area=0*units::mT/units::m*units::ms,
        double gradient_tolerance=1e-5);
    
    Regular(
        Species const & species_a, Quantity const & R1_b_or_T1_b,
        Vector3<Real> const & M0_a, Vector3<Real> const & M0_b,
        Quantity const & k_a,
        unsigned int initial_size=100, 
        Quantity const & unit_gradient_area=0*units::mT/units::m*units::ms,
        double gradient_tolerance=1e-5);
    
    Regular(Regular const &) = default;
    Regular(Regular &&) = default;
    Regular & operator=(Regular const &) = default;
    Regular & operator=(Regular &&) = default;
    virtual ~Regular() = default;
    
    /// @brief Return the number of states of the model.
    virtual std::size_t size() const;
    
    /// @brief Return the orders or the models.
    std::vector<Order> orders() const;
    
    using Base::state;
    
    /// @brief Return a given state of the model.
    ArrayC state(Order const & order) const;
    
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

    /**
     * @brief Simulate diffusion during given duration with given gradient
     * amplitude.
     */
    void diffusion(Quantity const & duration, Quantity const & gradient);
    
    void bulk_motion(Quantity const & duration, Quantity const & gradient);
    
    Quantity const & unit_gradient_area() const;
    double gradient_tolerance() const;
    
private:
    std::size_t _states_count;
    
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
