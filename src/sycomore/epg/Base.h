#ifndef _d635f223_7e0b_4f80_ab43_c03b499864f2
#define _d635f223_7e0b_4f80_ab43_c03b499864f2

#include <vector>

#include "sycomore/epg/Model.h"
#include "sycomore/magnetization.h"
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
    using State = std::vector<Complex>;
    using States = std::vector<Complex>;
    
    Quantity delta_omega=0*units::Hz;
    
    Real threshold=0;
    
    Base(
        Species const & species, Magnetization const & initial_magnetization,
        unsigned int initial_size);
    
    Base(
        Species const & species_a, Species const & species_b,
        Magnetization const & M0_a, Magnetization const & M0_b,
        Quantity const & k_a, Quantity const & delta_b,
        unsigned int initial_size);
    
    Base(
        Species const & species_a, Quantity const & R1_b_or_T1_b,
        Magnetization const & M0_a, Magnetization const & M0_b,
        Quantity const & k_a,
        unsigned int initial_size);
    
    Base(Base const & other) = default;
    Base(Base && other) = default;
    Base & operator=(Base const & other) = default;
    Base & operator=(Base && other) = default;
    virtual ~Base() = default;
    
    Species const & get_species(std::size_t pool=0) const;
    void set_species(std::size_t pool, Species const & species);
    void set_species(Species const & species);
    
    /// @brief Return the number of states of a pool.
    virtual std::size_t size() const = 0;
    
    /**
     * @brief Return a given state of the model, as a concatenation of
     * (F, F*, Z) for each pool.
     */
    State state(std::size_t order) const;
    
    /**
     * @brief Return all states in the model, as a concatenation of
     * (F, F*, Z) for each order and each pool.
     */
    States states() const;
    
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
};

}

}

#endif // _d635f223_7e0b_4f80_ab43_c03b499864f2
