#ifndef _8db2389d_b425_4fa0_8897_04a4ff117e15
#define _8db2389d_b425_4fa0_8897_04a4ff117e15

#include <xtensor/xtensor.hpp>

#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"
#include "sycomore/isochromat/Operator.h"

namespace sycomore
{

namespace isochromat
{

/// @brief Isochromat simulator
class SYCOMORE_API Model
{
public:
    /// @brief Create a spatially constant model
    Model(
        Quantity const & T1, Quantity const & T2, TensorR<1> const & M0, 
        TensorQ<2> const & positions,
        Quantity const & delta_omega=0*sycomore::units::Hz);
    
    /// @brief Create a spatially-varying model
    Model(
        TensorQ<1> const & T1, TensorQ<1> const & T2, TensorR<2> const & M0,
        TensorQ<2> const & positions,
        TensorQ<1> const & delta_omega={});
    
    /// @brief Create a spatially constant RF pulse operator
    Operator build_pulse(Quantity const & angle, Quantity const & phase) const;
    
    /// @brief Create a spatially-varying RF pulse operator
    Operator build_pulse(
        TensorQ<1> const & angle, TensorQ<1> const & phase) const;
    
    /// @brief Create a spatially constant time interval operator
    Operator build_time_interval(
        Quantity const & duration, Quantity const & delta_omega,
        TensorQ<1> const & gradient) const;
    
    /// @brief Create a spatially-varying time interval operator
    Operator build_time_interval(
        Quantity const & duration, TensorQ<1> const & delta_omega,
        TensorQ<2> const & gradient) const;
    
    /// @brief Create a relaxation operator
    Operator build_relaxation(Quantity const & duration) const;
    
    /// @brief Create a spatially constant phase accumulation operator
    Operator build_phase_accumulation(Quantity const & angle) const;
    
    /// @brief Create a spatially-varying phase accumulation operator
    Operator build_phase_accumulation(TensorQ<1> const & angle) const;
    
    /// @brief Apply an operator to the magnetization
    void apply(Operator const & operator_);
    
    /// @brief Return the T1 field
    TensorQ<1> T1() const;
    
    /// @brief Return the T2 field
    TensorQ<1> T2() const;
    
    /// @brief Return the M0 field
    TensorR<1> const & M0() const;
    
    /// @brief Return the off-resonance field.
    TensorR<1> const & delta_omega() const;
    
    /// @brief Return the magnetization field
    TensorR<2> magnetization() const;
    
    /// @brief Return the positions
    TensorQ<2> positions() const;
    
private:
    TensorR<1> _T1;
    TensorR<1> _T2;
    TensorR<1> _M0;
    TensorR<1> _delta_omega;
    
    TensorR<2> _magnetization;
    TensorR<2> _positions;
};

}

}

#endif // _8db2389d_b425_4fa0_8897_04a4ff117e15
