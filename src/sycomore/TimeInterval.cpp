#include "TimeInterval.h"

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

TimeInterval
::TimeInterval(Quantity const & duration, Quantity const & gradient_moment)
{
    this->set_duration(duration);
    this->set_gradient_moment({gradient_moment, gradient_moment, gradient_moment});
}

TimeInterval
::TimeInterval(
    Quantity const & duration, Array<Quantity> const & gradient_moment)
{
    this->set_duration(duration);
    this->set_gradient_moment(gradient_moment);
}

Quantity const &
TimeInterval
::get_duration() const
{
    return this->_duration;
}

void
TimeInterval
::set_duration(Quantity const & q)
{
    if(q.dimensions == Time)
    {
        this->_duration = q;
    }
    else
    {
        std::ostringstream message;
        message << "Invalid duration dimensions: " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

Array<Quantity> const &
TimeInterval
::get_gradient_moment() const
{
    return this->_gradient_moment;
}

void
TimeInterval
::set_gradient_moment(Array<Quantity> const & a)
{
    for(auto && q:a)
    {
        if(q.dimensions != GradientMoment)
        {
            std::ostringstream message;
            message << "Invalid gradient moment dimensions: " << q.dimensions;
            throw std::runtime_error(message.str());
        }
    }

    this->_gradient_moment = a;
}

bool
TimeInterval
::operator==(TimeInterval const & other) const
{
    return (
        this->_duration == other._duration
        && this->_gradient_moment == other._gradient_moment);
}

bool
TimeInterval
::operator!=(TimeInterval const & other) const
{
    return !(*this == other);
}

}
