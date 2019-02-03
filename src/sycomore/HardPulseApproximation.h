#ifndef _ea452652_2d33_45ed_9feb_960f1d861041
#define _ea452652_2d33_45ed_9feb_960f1d861041

#include <functional>
#include <string>
#include <vector>

#include "sycomore/Pulse.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

namespace sycomore
{

class HardPulseApproximation
{
public:
    using Envelope = std::function<Real(Quantity const &)>;

    HardPulseApproximation(
        Pulse const & model, std::vector<Quantity> const & support,
        Envelope const & envelope, std::string const & name);

    /**
     * @brief Create a hard pulse approximation with a slice-selection gradient
     * during the time interval.
     */
    HardPulseApproximation(
        Pulse const & model, std::vector<Quantity> const & support,
        Envelope const & envelope, Quantity const & bandwidth,
        Quantity const & slice_thickness, std::string const & name);

    std::vector<Pulse> const & get_pulses() const;
    TimeInterval const & get_time_interval() const;
    std::string const & get_name() const;

    Array<Real> get_gradient_moment() const;

    void set_phase(Quantity const & phase);
    void set_phase(Real const & phase);

private:
    std::vector<Pulse> _pulses;
    TimeInterval _time_interval;
    std::string _name;
};

HardPulseApproximation::Envelope sinc_envelope(Quantity const & t0);

}

#endif // _ea452652_2d33_45ed_9feb_960f1d861041