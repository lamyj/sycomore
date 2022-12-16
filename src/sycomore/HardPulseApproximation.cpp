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
    Envelope const & envelope)
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
    Quantity const & slice_thickness)
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

HardPulseApproximation::Envelope
apodized_sinc_envelope(Quantity const & t0_, unsigned int N, Real alpha)
{
    return [=](Quantity const t_) {
        auto const t0 = t0_.magnitude;
        auto const t = t_.magnitude;
        return 
            t0
            *((1-alpha)+alpha*std::cos(M_PI*t/(N*t0)))
            *std::sin(M_PI*t/t0)/(M_PI*t);
    };
}

HardPulseApproximation::Envelope
hanning_sinc_envelope(Quantity const & t0, unsigned int N)
{
    return apodized_sinc_envelope(t0, N, 0.5);
}

HardPulseApproximation::Envelope
hamming_sinc_envelope(Quantity const & t0, unsigned int N)
{
    return apodized_sinc_envelope(t0, N, 0.46);
}

HardPulseApproximation::Envelope sinc_envelope(Quantity const & t0)
{
    return [&](Quantity const x) {
        double const x_scaled = M_PI*x/t0;
        return x_scaled==0?1:std::sin(x_scaled)/(x_scaled);
    };
}

}
