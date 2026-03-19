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

// Gradient dephasing is in rad/m, i.e. [ L^-1 ]
Dimensions const GradientDephasing{-1, 0, 0, 0, 0, 0, 0};

TensorQ<1> linspace(Quantity const & min, Quantity const & max, std::size_t size)
{
    TensorQ<1> result(TensorQ<1>::shape_type{size});
    result.dimensions = min.dimensions;
    auto const delta = (max-min)/(size-1);
    for(std::size_t i=0; i<size; ++i)
    {
        result.unchecked(i) = min+i*delta;
    }
    return result;
}

TensorQ<1> linspace(Quantity const & span, std::size_t size)
{
    return linspace(-span/2., span/2., size);
}

// From CODATA 2018: https://www.physics.nist.gov/cgi-bin/cuu/Value?gammap
// rad*MHz/T simplifies to [ M^-1 * T * I] since T is [ M T^-2 I^-1 ]
Quantity const gamma{2.6752218744*1e8, {0, -1, 1, 1, 0, 0, 0}};
// https://www.physics.nist.gov/cgi-bin/cuu/Value?gammapbar
Quantity const gamma_bar{42.577478518*1e6, {0, -1, 1, 1, 0, 0, 0}};

}
