#ifndef _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
#define _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74

#include <vector>

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

struct TimeInterval
{
    /// @brief Interval duration in seconds
    Real duration;

    /// @brief Gradient moment on x,y,z axes
    std::vector<Real> gradient_moment;

    // TODO: gradient shape

    TimeInterval(Real duration, Real gradient_moment=0);
    TimeInterval(units::Time duration, Real gradient_moment=0);
    TimeInterval(Real duration, std::vector<Real> gradient_moment);
    TimeInterval(units::Time duration, std::vector<Real> gradient_moment);
};

}

#endif // _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
