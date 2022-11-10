#ifndef _8db2389d_b425_4fa0_8897_04a4ff117e15
#define _8db2389d_b425_4fa0_8897_04a4ff117e15

#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

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
        Real T1, Real T2,
        xt::xtensor<Real, 1> const & M0, xt::xtensor<Real, 2> const & positions);
    Model(
        xt::xtensor<Real, 1> const & T1, xt::xtensor<Real, 1> const & T2,
        xt::xtensor<Real, 2> const & M0, xt::xtensor<Real, 2> const & positions);
    
    Operator build_pulse(Real angle, Real phase) const;
    Operator build_pulse(
        xt::xtensor<Real, 1> const & angle,
        xt::xtensor<Real, 1> const & phase) const;
    
    Operator build_time_interval(
        Real duration, Real delta_omega,
        xt::xtensor<Real, 1> const & gradient) const;
    Operator build_time_interval(
        Real duration, xt::xtensor<Real, 1> const & delta_omega,
        xt::xtensor<Real, 2> const & gradient) const;
    
    Operator build_relaxation(Real duration) const;
    
    Operator build_phase_accumulation(Real angle) const;
    Operator build_phase_accumulation(xt::xtensor<Real, 1> const & angle) const;
    
    xt::xtensor<Real, 1> const & T1() const;
    xt::xtensor<Real, 1> const & T2() const;
    xt::xtensor<Real, 1> const & M0() const;
    xt::xtensor<Real, 2> const & magnetization() const;
    xt::xtensor<Real, 2> const & positions() const;
    
private:
    xt::xtensor<Real, 1> _T1;
    xt::xtensor<Real, 1> _T2;
    xt::xtensor<Real, 1> _M0;
    xt::xtensor<Real, 2> _magnetization;
    xt::xtensor<Real, 2> _positions;
};

}

}

#endif // _8db2389d_b425_4fa0_8897_04a4ff117e15
