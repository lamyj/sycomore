#include "Pulse.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

Pulse
::Pulse(Quantity const & angle, Quantity const & phase)
{
    this->set_angle(angle);
    this->set_phase(phase);
}

Quantity const &
Pulse
::get_angle() const
{
    return this->_angle;
}

void
Pulse
::set_angle(Quantity const & q)
{
    if(q.dimensions == Angle)
    {
        this->_angle = q;
    }
    else
    {
        std::ostringstream message;
        message << "Invalid angle dimensions: " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

Quantity const &
Pulse
::get_phase() const
{
    return this->_phase;
}

void
Pulse
::set_phase(Quantity const & q)
{
    if(q.dimensions == Angle)
    {
        this->_phase = q;
    }
    else
    {
        std::ostringstream message;
        message << "Invalid phase dimensions: " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

}
