#ifndef _a5b5eb59_d6dc_4067_b736_f03a01085d12
#define _a5b5eb59_d6dc_4067_b736_f03a01085d12

#include <cstdint>
#include <complex>
#include <vector>

#include <xtensor/xfixed.hpp>

#include "sycomore/Array.h"
#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

using Complex = std::complex<Real>;
using Magnetization = xt::xtensor_fixed<Real, xt::xshape<3>>;

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

template<typename T>
T round(T const & x, T const & r)
{
    return std::round(x/r)*r;
}

/// @brief Gyromagnetic ratio of 1H in rad/s/T
SYCOMORE_API extern Quantity const gamma;

/// @brief Gyromagnetic ratio of 1H in Hz/T
SYCOMORE_API extern Quantity const gamma_bar;

}

#endif // _a5b5eb59_d6dc_4067_b736_f03a01085d12
