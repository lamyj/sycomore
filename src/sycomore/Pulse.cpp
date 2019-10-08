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
        this->_angle= q;
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

Pulse::RotationMatrix
Pulse
::rotation_matrix() const
{
    using namespace units;

    Real const c_alpha = std::cos(this->_angle.convert_to(rad));
    Real const s_alpha = std::sin(this->_angle.convert_to(rad));
    Complex const e_i_phi{
        std::cos(this->_phase.convert_to(rad)),
        std::sin(this->_phase.convert_to(rad))};

    RotationMatrix m({0,0}, {3,3}, 0);
    m[{0, 0}] = 0.5 * (1. + c_alpha);
    m[{0, 1}] = Complex{0, -s_alpha} * e_i_phi / std::sqrt(2.);
    m[{0, 2}] = 0.5 * (1. - c_alpha) * std::pow(e_i_phi, 2);
    m[{1, 0}] = -std::conj(m[{0, 1}]);
    m[{1, 1}] = c_alpha;
    m[{1, 2}] = -m[{0, 1}];
    m[{2, 0}] = std::conj(m[{0, 2}]);
    m[{2, 1}] = -m[{1, 0}];
    m[{2, 2}] = m[{0, 0}];

    return m;
}

}
