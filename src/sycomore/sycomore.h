#ifndef _a5b5eb59_d6dc_4067_b736_f03a01085d12
#define _a5b5eb59_d6dc_4067_b736_f03a01085d12

#include <cstdint>
#include <complex>
#include <vector>

#if __has_include(<xtensor/xtensor.hpp>)
#include <xtensor/xfixed.hpp>
#else
#include <xtensor/containers/xfixed.hpp>
#endif

#include "sycomore/Array.h"
#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"
#include "sycomore/QuantityTensor.h"
#include "sycomore/sycomore_api.h"

/// @namespace sycomore @brief MRI Simulation Toolkit
namespace sycomore
{

/// @addtogroup Misc
/// @{

/// @brief Diffusion coefficient
SYCOMORE_API extern Dimensions const Diffusion;

/// @brief Gradient dephasing, as \f$\int \gamma G(t) dt\f$
SYCOMORE_API extern Dimensions const GradientDephasing;

/// @brief Generate evenly-spaced samples
TensorQ<1> linspace(Quantity const & min, Quantity const & max, std::size_t size);

/// @brief Generate evenly-spaced samples
TensorQ<1> linspace(Quantity const & span, std::size_t size);

/// @brief Round to x nearest multiple of r
template<typename T>
T round(T const & x, T const & r)
{
    return std::round(x/r)*r;
}

/// @brief Gyromagnetic ratio of 1H in rad/s/T
SYCOMORE_API extern Quantity const gamma;

/// @brief Gyromagnetic ratio of 1H in Hz/T
SYCOMORE_API extern Quantity const gamma_bar;

/// @}

}

#endif // _a5b5eb59_d6dc_4067_b736_f03a01085d12
