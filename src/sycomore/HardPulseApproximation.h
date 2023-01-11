#ifndef _ea452652_2d33_45ed_9feb_960f1d861041
#define _ea452652_2d33_45ed_9feb_960f1d861041

#include <functional>
#include <string>
#include <vector>

#include "sycomore/Pulse.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

namespace sycomore
{

class SYCOMORE_API HardPulseApproximation
{
public:
    using Envelope = std::function<Quantity(Quantity const &)>;

    HardPulseApproximation(
        Pulse const & model, std::vector<Quantity> const & support,
        Envelope const & envelope);

    std::vector<Pulse> const & get_pulses() const;
    TimeInterval const & get_time_interval() const;

    Vector3Q get_gradient_moment() const;

    void set_phase(Quantity const & phase);

private:
    std::vector<Pulse> _pulses;
    TimeInterval _time_interval;
};

SYCOMORE_API HardPulseApproximation::Envelope apodized_sinc_envelope(
    Quantity const & t0, unsigned int N, Real alpha);
SYCOMORE_API HardPulseApproximation::Envelope hanning_sinc_envelope(
    Quantity const & t0, unsigned int N);
SYCOMORE_API HardPulseApproximation::Envelope hamming_sinc_envelope(
    Quantity const & t0, unsigned int N);
SYCOMORE_API HardPulseApproximation::Envelope sinc_envelope(Quantity const & t0);

}

#endif // _ea452652_2d33_45ed_9feb_960f1d861041
