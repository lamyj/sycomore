#include "TimeInterval.h"

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

TimeInterval
TimeInterval
::shortest(Quantity const & k, Quantity const & G_max)
{
    return TimeInterval::shortest({k, k, k}, G_max);
}

TimeInterval
TimeInterval
::shortest(Vector3Q const & k, Quantity const & G_max)
{
    Quantity max(xt::amax(xt::abs(k.magnitude))(), k.dimensions);
    
    if(max.dimensions == (units::rad/units::m).dimensions)
    {
        max /= gamma;
    }
    if(max.dimensions != (units::T/units::m*units::s).dimensions)
    {
        throw std::runtime_error("Invalid dimensions");
    }
    
    auto const min_time = max/G_max;
    return sycomore::TimeInterval(min_time, k);
}

TimeInterval
::TimeInterval()
: TimeInterval(0*sycomore::units::s, 0*sycomore::units::T/sycomore::units::m)
{
    // Nothing else.
}

TimeInterval
::TimeInterval(
    Quantity const & duration, Quantity const & gradient)
: TimeInterval(duration, {gradient, gradient, gradient})
{
    // Nothing else.
}

TimeInterval
::TimeInterval(Quantity const & duration, Vector3Q const & gradient)
{
    this->set_duration(duration);
    this->set_gradient(gradient);
}

Quantity const &
TimeInterval
::duration() const
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

void
TimeInterval
::set_gradient(Quantity const & q)
{
    this->set_gradient({q,q,q});
}

void
TimeInterval
::set_gradient(Vector3Q const & a)
{
    if(a.dimensions == (units::T/units::m).dimensions)
    {
        this->_gradient_amplitude = a;
    }
    else if(a.dimensions == (units::T/units::m*units::s).dimensions)
    {
        if(this->_duration == 0*units::s)
        {
            this->_gradient_amplitude.fill(0*units::T/units::m);
        }
        else
        {
            this->_gradient_amplitude = a/this->_duration;
        }
    }
    else if(a.dimensions == (units::rad/units::m).dimensions)
    {
        if(this->_duration == 0*units::s)
        {
            this->_gradient_amplitude.fill(0*units::T/units::m);
        }
        else
        {
            this->_gradient_amplitude = a/(this->_duration*sycomore::gamma);
        }
    }
    else
    {
        std::ostringstream message;
        message << "Invalid gradient specification: " << a.dimensions;
        throw std::runtime_error(message.str());
    }
}

Vector3Q const &
TimeInterval
::gradient_amplitude() const
{
    return this->_gradient_amplitude;
}

Vector3Q
TimeInterval
::gradient_area() const
{
    return this->_duration*this->_gradient_amplitude;
}

Vector3Q
TimeInterval
::gradient_dephasing() const
{
    return sycomore::gamma*this->_duration*this->_gradient_amplitude;
}

bool
TimeInterval
::operator==(TimeInterval const & other) const
{
    return (
        this->_duration == other._duration
        && this->_gradient_amplitude == other._gradient_amplitude);
}

bool
TimeInterval
::operator!=(TimeInterval const & other) const
{
    return !(*this == other);
}

}
