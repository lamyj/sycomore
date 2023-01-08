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
class SYCOMORE_API Base
{
public:
    Quantity delta_omega=0*units::Hz;
    
    Real threshold=0;
    
    Base(
        Species const & species, Vector3<Real> const & initial_magnetization,
        unsigned int initial_size);
    
    Base(
        Species const & species_a, Species const & species_b,
        Vector3<Real> const & M0_a, Vector3<Real> const & M0_b,
        Quantity const & k_a, Quantity const & delta_b,
        unsigned int initial_size);
    
    Base(
        Species const & species_a, Quantity const & R1_b_or_T1_b,
        Vector3<Real> const & M0_a, Vector3<Real> const & M0_b,
        Quantity const & k_a,
        unsigned int initial_size);
    
    Base(Base const & other) = default;
    Base(Base && other) = default;
    Base & operator=(Base const & other) = default;
    Base & operator=(Base && other) = default;
    virtual ~Base() = default;
    
    Model::Kind kind() const;
    std::size_t pools() const;
    
    Species const & get_species(std::size_t pool=0) const;
    void set_species(std::size_t pool, Species const & species);
    void set_species(Species const & species);
    
    Real const & get_M0(std::size_t pool=0) const;
    void set_M0(std::size_t pool, Real const & M0);
    void set_M0(Real const & M0);
    
    Quantity const & get_k(std::size_t pool) const;
    void set_k(std::size_t pool, Quantity const & k);
    
    Quantity const & get_delta_b() const;
    void set_delta_b(Quantity const & delta_b);
    
    /// @brief Return the number of states of a pool.
    virtual std::size_t size() const = 0;
    
    /**
     * @brief Return a given state of the model, as a concatenation of
     * (F, F*, Z) for each pool.
     */
    ArrayC state(std::size_t order) const;
    
    /**
     * @brief Return all states in the model, as a concatenation of
     * (F, F*, Z) for each order and each pool.
     */
    ArrayC states() const;
    
    /// @brief Return the elapsed time.
    Quantity elapsed() const;
    
    /// @brief Return the echo signal, i.e. F_0
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
    Model _model;
    Real _elapsed;
};

}

}

#endif // _d635f223_7e0b_4f80_ab43_c03b499864f2
