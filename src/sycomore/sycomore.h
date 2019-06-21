#ifndef _a5b5eb59_d6dc_4067_b736_f03a01085d12
#define _a5b5eb59_d6dc_4067_b736_f03a01085d12

#include <cstdint>
#include <complex>
#include <vector>

#include "sycomore/Array.h"
#include "sycomore/units.h"

namespace sycomore
{

using Real = double;
using Complex = std::complex<Real>;

using Index = Array<int32_t>;
using Shape = Array<uint32_t>;
using Stride = Array<uint32_t>;

using Point = Array<Quantity>;

auto const Diffusion = std::pow(Length, 2)/Time;

// NOTE: gradient moment is gamma * integral(g(t)*dt), with gamma in
// [rad * T^-1 * s^-1] and g(t) in [T * m^-1] (thus integral(g(t)*dt)
// in [T * m^-1 * s]). Gradient moment is then in [rad * m^-1], i.e. accumulated
// phase per length.
auto const GradientMoment = Angle/Length;

template<typename T>
std::vector<T> linspace(T min, T max, size_t size)
{
    std::vector<T> result;
    result.reserve(size);
    auto const delta = (max-min)/(size-1);
    for(size_t i=0; i<size; ++i)
    {
        result.push_back(min+i*delta);
    }
    return result;
}

template<typename T>
std::vector<T> linspace(T span, size_t size)
{
    return linspace(-span/2., span/2., size);
}

/// @brief Gyromagnetic ratio of 1H
Quantity const gamma=2*M_PI*units::rad * 42.57747892*units::MHz/units::T;

}

#endif // _a5b5eb59_d6dc_4067_b736_f03a01085d12
