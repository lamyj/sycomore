#ifndef _8db2389d_b425_4fa0_8897_04a4ff117e15
#define _8db2389d_b425_4fa0_8897_04a4ff117e15

#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

#include "sycomore/sycomore.h"

namespace sycomore
{

namespace isochromat
{

using Position = xt::xtensor_fixed<Real, xt::xshape<3>>;
using Positions = xt::xtensor<Position::value_type, 2>;
using Magnetization = xt::xtensor_fixed<Real, xt::xshape<3>>;
using Magnetizations = xt::xtensor<Magnetization::value_type, 2>;

using Operator = xt::xtensor<Real, 3>;

class SYCOMORE_API Model
{
public:
    using Array = xt::xtensor<Real, 1>;
    
    Model(Real T1, Real T2, Magnetization M0, Positions const & positions);
    Model(
        Array const & T1, Array const & T2, Magnetizations const & M0,
        Positions const & positions);
    
    Operator build_pulse(Real angle, Real phase) const;
    Operator build_pulse(Array const & angle, Array const & phase) const;
    void build_pulse(Real angle, Real phase, Operator & op) const;
    void build_pulse(
        Array const & angle, Array const & phase, Operator & op) const;
    
    Operator build_relaxation(Real duration) const;
    void build_relaxation(Real duration, Operator & op) const;

private:
    using Matrix = xt::xtensor_fixed<Operator::value_type, xt::xshape<4, 4>>;
    
    Array _T1;
    Array _T2;
    Array _M0;
    Magnetizations _magnetization;
    Positions _positions;
    
    void _combine(Operator const & left, Operator & right) const;
};

}

}

#endif // _8db2389d_b425_4fa0_8897_04a4ff117e15
