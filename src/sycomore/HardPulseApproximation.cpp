#include "HardPulseApproximation.h"

#include <functional>
#include <string>
#include <vector>

#include "sycomore/Pulse.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

namespace sycomore
{

HardPulseApproximation
::HardPulseApproximation(
    Pulse const & model, std::vector<Quantity> const & support,
    Envelope const & envelope, std::string const & name)
: _name(name)
{
    std::vector<Quantity> angles(support.size());
    std::transform(support.begin(), support.end(), angles.begin(), envelope);
    auto const sum = std::accumulate(angles.begin(), angles.end(), 0.*units::rad);
    std::transform(
        angles.begin(), angles.end(), angles.begin(),
        [&](Quantity x) { return x*model.get_angle() / sum; });

    for(auto && angle: angles)
    {
        this->_pulses.emplace_back(angle, model.get_phase());
    }

    auto const pulse_duration = support.back()-support.front();
    this->_time_interval.set_duration(pulse_duration/(support.size()-1));
}

HardPulseApproximation
::HardPulseApproximation(
    Pulse const & model, std::vector<Quantity> const & support,
    Envelope const & envelope, Quantity const & bandwidth,
    Quantity const & slice_thickness, std::string const & name)
: HardPulseApproximation(model, support, envelope, name)
{
    // From Handbook, eq. 8.53, which gives gradient amplitude. Assuming a
    // constant gradient amplitude, we get the moment by multiplying by pulse
    // duration
    using namespace units;
    auto const pulse_duration = support.back()-support.front();
    auto const total_moment =
        2*M_PI*rad * bandwidth / slice_thickness * pulse_duration;
    this->_time_interval.set_gradient_moment(
        {0.*rad/m, 0.*rad/m, total_moment/(support.size()-1)});
}

std::vector<Pulse> const &
HardPulseApproximation
::get_pulses() const
{
    return this->_pulses;
}

TimeInterval const &
HardPulseApproximation
::get_time_interval() const
{
    return this->_time_interval;
}

std::string const &
HardPulseApproximation
::get_name() const
{
    return this->_name;
}

Array<Quantity>
HardPulseApproximation
::get_gradient_moment() const
{
    return this->_time_interval.get_gradient_moment() * (this->_pulses.size()-1);
}

void
HardPulseApproximation
::set_phase(Quantity const & phase)
{
    for(auto && pulse: this->_pulses)
    {
        pulse.set_phase(phase);
    }
}

HardPulseApproximation::Envelope sinc_envelope(Quantity const & t0)
{
    return [&](Quantity const x) {
        auto const x_scaled = x.convert_to(units::s)/t0.convert_to(units::s);
        return units::rad*(x_scaled==0?1:std::sin(x_scaled*M_PI)/(x_scaled*M_PI));
    };
}

}
