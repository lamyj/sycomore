#include "TimeInterval.h"

#include <valarray>

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
::TimeInterval(units::Time duration, Real gradient_moment)
: TimeInterval(duration.convert_to(units::s), gradient_moment)
{
    // Nothing else.
}

TimeInterval
::TimeInterval(Real duration, std::vector<Real> gradient_moment)
: duration(duration), gradient_moment(gradient_moment)
{
    // Nothing else.
}

TimeInterval
::TimeInterval(units::Time duration, std::vector<Real> gradient_moment)
: TimeInterval(duration.convert_to(units::s), gradient_moment)
{
    // Nothing else.
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
