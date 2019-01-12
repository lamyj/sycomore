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

}

#endif // _a5b5eb59_d6dc_4067_b736_f03a01085d12
