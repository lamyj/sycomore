#ifndef _d635f223_7e0b_4f80_ab43_c03b499864f2
#define _d635f223_7e0b_4f80_ab43_c03b499864f2

#include "sycomore/Array.h"
#include "sycomore/epg/Model.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

/// @brief Base class for all EPG models.
class Base
{
public:
    /// @brief Frequency offset of the simulator
    Quantity delta_omega=0*units::Hz;
    
    /// @brief Threshold used to cull states with low population
    Real threshold=0;
    
    /// @brief Create a single-pool model
    Base(
        Species const & species, Vector3R const & initial_magnetization,
        unsigned int initial_size);
    
    /// @brief Create an exchange model
    Base(
        Species const & species_a, Species const & species_b,
        Vector3R const & M0_a, Vector3R const & M0_b,
        Quantity const & k_a, Quantity const & delta_b,
        unsigned int initial_size);
    
    /// @brief Create an MT model
    Base(
        Species const & species_a, Quantity const & R1_b_or_T1_b,
        Vector3R const & M0_a, Vector3R const & M0_b,
        Quantity const & k_a,
        unsigned int initial_size);
    
    /// @brief Default copy constructor
    Base(Base const & other) = default;
    /// @brief Default move constructor
    Base(Base && other) = default;
    /// @brief Default copy assignment
    Base & operator=(Base const & other) = default;
    /// @brief Default move assignment
    Base & operator=(Base && other) = default;
    /// @brief Default destructor
    virtual ~Base() = default;
    
    /// @brief Return the kind of the model, set at creation
    Model::Kind kind() const;
    
    /// @brief Return the number of pools of the model, set at creation
    std::size_t pools() const;
    
    /// @brief Return the species of one of the pools
    Species const & species(std::size_t pool=0) const;
    
    /// @brief Set the species of one of the pools
    void set_species(std::size_t pool, Species const & species);
    
    /// @brief Set the species of the first pool
    void set_species(Species const & species);
    
    /// @brief Return the equilibrium magnetization of one of the pools
    Real const & M0(std::size_t pool=0) const;
    
    /// @brief Set the equilibrium magnetization of one of the pools
    void set_M0(std::size_t pool, Real const & M0);
    
    /// @brief Set the equilibrium magnetization of the first pool
    void set_M0(Real const & M0);
    
    /// @brief Return the exchange constant of one of the pools
    Quantity const & k(std::size_t pool) const;
    
    /// @brief Set the exchange constant of one of the pools
    void set_k(std::size_t pool, Quantity const & k);
    
    /// @brief Return the frequency offset of the second pool
    Quantity const & delta_b() const;
    
    /// @brief Set the frequency offset of the second pool
    void set_delta_b(Quantity const & delta_b);
    
    /// @brief Return the number of states of a pool.
    virtual std::size_t size() const = 0;
    
    /**
     * @brief Return a given state of the model, as a concatenation of
     * (F, F*, Z) for each pool.
     */
    ArrayC state(std::size_t order) const;
    
    /**
     * @brief Return all states in the model, as \f$\tilde{F}\f$,
     * \f$\tilde{F}^*\f$, and \f$\tilde{Z}\f$ for each order and each pool.
     */
    ArrayC states() const;
    
    /// @brief Return the elapsed time.
    Quantity elapsed() const;
    
    /// @brief Return the echo signal, i.e. \f$F_0\f$
    Complex const & echo(std::size_t pool=0) const;
    
    /// @brief Apply an RF hard pulse to a single-pool model.
    void apply_pulse(Quantity const & angle, Quantity const & phase=0*units::rad);
    
    /// @brief Apply an RF hard pulse to an two-pools exchange model.
    void apply_pulse(
        Quantity const & angle_a, Quantity const & phase_a,
        Quantity const & angle_b, Quantity const & phase_b);
    
    /// @brief Apply an RF pulse to an two-pools magnetization transfer model.
    void apply_pulse(
        Quantity const & angle_a, Quantity const & phase_a, Real saturation);
    
    /// @brief Simulate the relaxation during given duration.
    void relaxation(Quantity const & duration);
    
    /**
     * @brief Simulate field- and species-related off-resonance effects during 
     * given duration with given frequency offset.
     */
    void off_resonance(Quantity const & duration);
    
protected:
    /// @brief EPG model
    Model _model;
    
    /// @brief Elapsed time, in s
    Real _elapsed;
};

}

}

#endif // _d635f223_7e0b_4f80_ab43_c03b499864f2
