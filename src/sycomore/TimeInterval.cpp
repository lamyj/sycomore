#include "TimeInterval.h"

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

TimeInterval
::TimeInterval(
    Quantity const & duration, Quantity const & gradient)
: TimeInterval(duration, {gradient, gradient, gradient})
{
    // Nothing else.
}

TimeInterval
::TimeInterval(Quantity const & duration, Array<Quantity> const & gradient)
{
    this->set_duration(duration);
    this->set_gradient(gradient);
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

void
TimeInterval
::set_gradient(Quantity const & q)
{
    this->set_gradient({q,q,q});
}

void
TimeInterval
::set_gradient(Array<Quantity> const & a)
{
    if(a.empty())
    {
        throw std::runtime_error("Cannot set gradient from empty array");
    }
    auto const & q = a[0];
    
    if(q.dimensions == (units::T/units::m).dimensions)
    {
        this->set_gradient_amplitude(a);
    }
    else if(q.dimensions == (units::T/units::m*units::s).dimensions)
    {
        this->set_gradient_area(a);
    }
    else if(q.dimensions == (units::rad/units::m).dimensions)
    {
        this->set_gradient_dephasing(a);
    }
    else
    {
        std::ostringstream message;
        message << "Invalid gradient specification: " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

Array<Quantity>
TimeInterval
::get_gradient_moment() const
{
    return sycomore::gamma*this->_duration*this->_gradient_amplitude;
}

void
TimeInterval
::set_gradient_moment(Quantity const & q)
{
    this->set_gradient_moment({q,q,q});
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
    
    if(this->_duration == 0*units::s)
    {
        this->set_gradient_amplitude(0*units::T/units::m);
    }
    else
    {
        this->set_gradient_amplitude(a/(this->_duration*sycomore::gamma));
    }
}

Array<Quantity> const &
TimeInterval
::get_gradient_amplitude() const
{
    return this->_gradient_amplitude;
}

void
TimeInterval
::set_gradient_amplitude(Quantity const & q)
{
    this->set_gradient_amplitude({q,q,q});
}

void
TimeInterval
::set_gradient_amplitude(Array<Quantity> const & a)
{
    for(auto && q:a)
    {
        if(q.dimensions != (units::T/units::m).dimensions)
        {
            std::ostringstream message;
            message << "Invalid gradient amplitude dimensions: " << q.dimensions;
            throw std::runtime_error(message.str());
        }
    }

    this->_gradient_amplitude = a;
}

Array<Quantity>
TimeInterval
::get_gradient_area() const
{
    return this->_duration*this->_gradient_amplitude;
}

void
TimeInterval
::set_gradient_area(Quantity const & q)
{
    this->set_gradient_area({q,q,q});
}

void
TimeInterval
::set_gradient_area(Array<Quantity> const & a)
{
    for(auto && q:a)
    {
        if(q.dimensions != (units::T/units::m*units::s).dimensions)
        {
            std::ostringstream message;
            message << "Invalid gradient area dimensions: " << q.dimensions;
            throw std::runtime_error(message.str());
        }
    }

    if(this->_duration == 0*units::s)
    {
        this->set_gradient_amplitude(0*units::T/units::m);
    }
    else
    {
        this->set_gradient_amplitude(a/this->_duration);
    }
}

Array<Quantity>
TimeInterval
::get_gradient_dephasing() const
{
    return this->get_gradient_moment();
}

void
TimeInterval
::set_gradient_dephasing(Quantity const & q)
{
    this->set_gradient_dephasing({q,q,q});
}

void
TimeInterval
::set_gradient_dephasing(Array<Quantity> const & a)
{
    this->set_gradient_moment(a);
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
