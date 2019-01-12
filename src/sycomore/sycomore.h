#ifndef _a5b5eb59_d6dc_4067_b736_f03a01085d12
#define _a5b5eb59_d6dc_4067_b736_f03a01085d12

#include <complex>

#include "sycomore/Array.h"
#include "sycomore/units.h"

namespace sycomore
{

using Real = double;
using Complex = std::complex<Real>;

using Index = Array<int>;
using Shape = Array<unsigned int>;
using Stride = Array<unsigned int>;

using Diffusion = units::div<units::pow<units::Length, 2>, units::Time>;

// NOTE: gradient moment is gamma * integral(g(t)*dt), with gamma in
// [T^-1 * s^-1] and g(t) in [T * m^-1] (thus integral(g(t)*dt)
// in [T * m^-1 * s]). Gradient moment is then in [m^-1]
using GradientMoment = units::pow<units::Length, -1>;

}

#endif // _a5b5eb59_d6dc_4067_b736_f03a01085d12
