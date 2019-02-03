#ifndef _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
#define _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

struct TimeInterval
{
    /// @brief Interval duration in seconds
    Real duration;

    /// @brief Gradient moment on x,y,z axes, in 1/m.
    Array<Real> gradient_moment;

    // TODO: gradient shape

    TimeInterval(Real duration=0, Real gradient_moment=0);
    TimeInterval(
        Quantity const & duration, Quantity const & gradient_moment={0,GradientMoment});
    TimeInterval(Real duration, Array<Real> gradient_moment);
    TimeInterval(Quantity const & duration, Array<Quantity> const & gradient_moment);

    bool operator==(TimeInterval const & other) const;
    bool operator!=(TimeInterval const & other) const;
};

}

#endif // _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
