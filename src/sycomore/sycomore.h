#ifndef _a5b5eb59_d6dc_4067_b736_f03a01085d12
#define _a5b5eb59_d6dc_4067_b736_f03a01085d12

#include <cstdint>
#include <complex>
#include <vector>

#include "sycomore/Array.h"
#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

using Real = double;
using Complex = std::complex<Real>;

using Index = Array<int32_t>;
using Shape = Array<uint32_t>;
using Stride = Array<uint32_t>;

using Point = Array<Quantity>;

SYCOMORE_API extern Dimensions const Diffusion;

// NOTE: gradient moment is gamma * integral(g(t)*dt), with gamma in
// [rad * T^-1 * s^-1] and g(t) in [T * m^-1] (thus integral(g(t)*dt)
// in [T * m^-1 * s]). Gradient moment is then in [rad * m^-1], i.e. accumulated
// phase per length.
SYCOMORE_API extern Dimensions const GradientMoment;

template<typename T>
std::vector<T> linspace(T min, T max, std::size_t size)
{
    std::vector<T> result;
    result.reserve(size);
    auto const delta = (max-min)/(size-1);
    for(std::size_t i=0; i<size; ++i)
    {
        result.push_back(min+i*delta);
    }
    return result;
}

template<typename T>
std::vector<T> linspace(T span, std::size_t size)
{
    return linspace(-span/2., span/2., size);
}

/// @brief Gyromagnetic ratio of 1H
SYCOMORE_API extern Quantity const gamma;

}

#endif // _a5b5eb59_d6dc_4067_b736_f03a01085d12
