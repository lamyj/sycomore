#include "HardPulseApproximation.h"

#include <functional>
#include <string>
#include <vector>

#include "sycomore/Pulse.h"
#include "sycomore/Quantity.h"
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
        [&](Quantity x) { return x*model.angle() / sum; });

    for(auto && angle: angles)
    {
        this->_pulses.emplace_back(angle, model.phase());
    }

    auto const pulse_duration = support.back()-support.front();
    this->_duration = pulse_duration/(support.size()-1);
}

std::vector<Pulse> const &
HardPulseApproximation
::pulses() const
{
    return this->_pulses;
}

Quantity const &
HardPulseApproximation
::duration() const
{
    return this->_duration;
}

Quantity const &
HardPulseApproximation
::phase() const
{
    return this->_pulses.front().phase();
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
