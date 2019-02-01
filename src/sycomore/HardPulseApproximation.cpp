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
    Pulse const & model, std::vector<units::Time> const & support,
    Envelope const & envelope, std::string const & name)
: _name(name)
{
    std::vector<Real> angles(support.size());
    std::transform(support.begin(), support.end(), angles.begin(), envelope);
    auto const sum = std::accumulate(angles.begin(), angles.end(), 0.);
    std::transform(
        angles.begin(), angles.end(), angles.begin(),
        [&](Real x) { return x*model.angle / sum; });

    for(auto && angle: angles)
    {
        this->_pulses.emplace_back(angle, model.phase);
    }

    auto const pulse_duration = support.back()-support.front();
    this->_time_interval.duration = pulse_duration.value/(support.size()-1);
}

HardPulseApproximation
::HardPulseApproximation(
    Pulse const & model, std::vector<units::Time> const & support,
    Envelope const & envelope, units::Frequency const & bandwidth,
    units::Length const & slice_thickness, std::string const & name)
: HardPulseApproximation(model, support, envelope, name)
{
    // From Handbook, eq. 8.53, which gives gradient amplitude. Assuming a
    // constant gradient amplitude, we get the moment by multiplying by pulse
    // duration
    auto const pulse_duration = support.back()-support.front();
    auto const total_moment =
        2*M_PI*bandwidth.value/slice_thickness.value * pulse_duration.value;
    this->_time_interval.gradient_moment = {
        0., 0., total_moment/(support.size()-1)};
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

Array<Real>
HardPulseApproximation
::get_gradient_moment() const
{
    return this->_time_interval.gradient_moment * (this->_pulses.size()-1);
}

void
HardPulseApproximation
::set_phase(units::Angle const & phase)
{
    this->set_phase(phase.convert_to(units::rad));
}

void
HardPulseApproximation
::set_phase(Real const & phase)
{
    for(auto && pulse: this->_pulses)
    {
        pulse.phase = phase;
    }
}

HardPulseApproximation::Envelope sinc_envelope(units::Time const & t0)
{
    return [&](units::Time const x) {
        auto const x_scaled = x.value/t0.value;
        return x_scaled==0?1:std::sin(x_scaled*M_PI)/(x_scaled*M_PI);
    };
}

}
