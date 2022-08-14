#ifndef _d635f223_7e0b_4f80_ab43_c03b499864f2
#define _d635f223_7e0b_4f80_ab43_c03b499864f2

#include <memory>

#include "sycomore/epg/pool_model.h"
#include "sycomore/epg/pool_storage.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

/// @brief Base class for all EPG models.
class Base
{
public:
    Quantity delta_omega=0*units::Hz;
    
    Real threshold=0;
    
    Base(
        Species const & species, Magnetization const & initial_magnetization,
        unsigned int initial_size);
    Base(Base const & other);
    Base(Base && other);
    Base & operator=(Base const & other);
    Base & operator=(Base && other);
    virtual ~Base() = default;
    
    Species const & get_species() const;
    void set_species(Species const & species);
    
    virtual std::size_t size() const = 0;
    
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
    
    /// @brief Simulate the relaxation during given duration.
    void relaxation(Quantity const & duration);
    
    /**
     * @brief Simulate field- and species-related off-resonance effects during 
     * given duration with given frequency offset.
     */
    void off_resonance(Quantity const & duration);
    
protected:
    std::shared_ptr<pool_storage::Base> _storage;
    std::shared_ptr<pool_model::Base> _model;
};

}

}

#endif // _d635f223_7e0b_4f80_ab43_c03b499864f2
