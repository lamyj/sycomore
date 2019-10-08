#include "sycomore.h"

#include <cmath>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"
#include "sycomore/units.h"

namespace sycomore
{

// WARNING 3.6.2-2: the order of initialization in different translation units 
// is undefined. Avoid using pre-defined dimensions or units.

// [ L^2 * T^-1 ]
Dimensions const Diffusion{2, 0, -1, 0, 0, 0, 0};

// Gradient moment is in rad/m, i.e. accumulated phase per length, i.e. [ L^-1 ]
Dimensions const GradientMoment{-1, 0, 0, 0, 0, 0, 0};

// rad*MHz/T simplifies to [ M^-1 * T * I] since T is [ M T^-2 I^-1 ]
Quantity const gamma{2*M_PI * 42.57747892*1e6, {0, -1, 1, 1, 0, 0, 0}};

}
