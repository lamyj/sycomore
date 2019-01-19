#include "TimeInterval.h"

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

TimeInterval
::TimeInterval(Real duration, Real gradient_moment)
: TimeInterval(duration, {gradient_moment, gradient_moment, gradient_moment})
{
    // Nothing else.
}

TimeInterval
::TimeInterval(units::Time duration, GradientMoment gradient_moment)
: TimeInterval(duration.convert_to(units::s),
    gradient_moment.convert_to(1/units::m))
{
    // Nothing else.
}

TimeInterval
::TimeInterval(Real duration, Array<Real> gradient_moment)
: duration(duration), gradient_moment(gradient_moment)
{
    // Nothing else.
}

TimeInterval
::TimeInterval(units::Time duration, Array<GradientMoment> gradient_moment_)
: duration(duration.convert_to(units::s))
{
    this->gradient_moment = Array<Real>(gradient_moment_.size());
    for(size_t d=0; d<this->gradient_moment.size(); ++d)
    {
        this->gradient_moment[d] = gradient_moment_[d].convert_to(1/units::m);
    }
}

bool
TimeInterval
::operator==(TimeInterval const & other) const
{
    return (
        this->duration == other.duration
        && this->gradient_moment == other.gradient_moment);
}

bool
TimeInterval
::operator!=(TimeInterval const & other) const
{
    return !(*this == other);
}

}
