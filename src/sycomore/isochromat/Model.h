#ifndef _8db2389d_b425_4fa0_8897_04a4ff117e15
#define _8db2389d_b425_4fa0_8897_04a4ff117e15

#include <xtensor/xtensor.hpp>

#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/isochromat/Operator.h"

namespace sycomore
{

namespace isochromat
{

class SYCOMORE_API Model
{
public:
    Model(
        Quantity const & T1, Quantity const & T2, TensorR<1> const & M0, 
        TensorQ<2> const & positions);
    Model(
        TensorQ<1> const & T1, TensorQ<1> const & T2, TensorR<2> const & M0,
        TensorQ<2> const & positions);
    
    Operator build_pulse(Quantity const & angle, Quantity const & phase) const;
    Operator build_pulse(
        TensorQ<1> const & angle, TensorQ<1> const & phase) const;
    
    Operator build_time_interval(
        Quantity const & duration, Quantity const & delta_omega,
        TensorQ<1> const & gradient) const;
    Operator build_time_interval(
        Quantity const & duration, TensorQ<1> const & delta_omega,
        TensorQ<2> const & gradient) const;
    
    Operator build_relaxation(Quantity const & duration) const;
    
    Operator build_phase_accumulation(Quantity const & angle) const;
    Operator build_phase_accumulation(TensorQ<1> const & angle) const;
    
    void apply(Operator const & operator_);
    
    TensorQ<1> T1() const;
    TensorQ<1> T2() const;
    TensorR<1> const & M0() const;
    TensorR<2> magnetization() const;
    TensorQ<2> positions() const;
    
private:
    TensorR<1> _T1;
    TensorR<1> _T2;
    TensorR<1> _M0;
    TensorR<2> _magnetization;
    TensorR<2> _positions;
};

}

}

#endif // _8db2389d_b425_4fa0_8897_04a4ff117e15
