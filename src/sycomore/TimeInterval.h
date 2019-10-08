#ifndef _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
#define _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

namespace sycomore
{

struct SYCOMORE_API TimeInterval
{
    // TODO: gradient shape

    TimeInterval(
        Quantity const & duration=0*units::s,
        Quantity const & gradient_moment=0*units::rad/units::m);
    TimeInterval(
        Quantity const & duration,
        Array<Quantity> const & gradient_moment);

    Quantity const & get_duration() const;
    void set_duration(Quantity const & q);

    Array<Quantity> const & get_gradient_moment() const;
    void set_gradient_moment(Array<Quantity> const & a);

    bool operator==(TimeInterval const & other) const;
    bool operator!=(TimeInterval const & other) const;

private:
    /// @brief Interval duration
    Quantity _duration;

    /// @brief Gradient moment on x,y,z axes.
    Array<Quantity> _gradient_moment;
};

}

#endif // _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
