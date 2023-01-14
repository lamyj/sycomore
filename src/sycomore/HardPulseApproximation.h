#ifndef _ea452652_2d33_45ed_9feb_960f1d861041
#define _ea452652_2d33_45ed_9feb_960f1d861041

#include <functional>
#include <string>
#include <vector>

#include "sycomore/Pulse.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

namespace sycomore
{

/// @brief Small tip angle approximation of a shaped pulse
class SYCOMORE_API HardPulseApproximation
{
public:
    /// @brief Normalized envelope of the pulse as a function of time
    using Envelope = std::function<Quantity(Quantity const &)>;
    
    /**
     * @brief Create a shaped pulse with a flip angle equivalent to given
     * hard pulse
     */
    HardPulseApproximation(
        Pulse const & model, std::vector<Quantity> const & support,
        Envelope const & envelope);
    
    /// @brief Return the hard pulses approximating the shaped pulse
    std::vector<Pulse> const & get_pulses() const;
    
    /// @brief Return the duration of a hard pulse
    Quantity const & get_duration() const;
    
    /// @brief Set the phase of the shaped pulse
    void set_phase(Quantity const & phase);

private:
    std::vector<Pulse> _pulses;
    Quantity _duration;
};

/// @brief Create an apodized sinc envelope
SYCOMORE_API HardPulseApproximation::Envelope apodized_sinc_envelope(
    Quantity const & t0, unsigned int N, Real alpha);
/// @brief Create an Hanning-apodized sinc envelope
SYCOMORE_API HardPulseApproximation::Envelope hanning_sinc_envelope(
    Quantity const & t0, unsigned int N);
/// @brief Create an Hamming-apodized sinc envelope
SYCOMORE_API HardPulseApproximation::Envelope hamming_sinc_envelope(
    Quantity const & t0, unsigned int N);
/// @brief Create a sinc envelope
SYCOMORE_API HardPulseApproximation::Envelope sinc_envelope(Quantity const & t0);

}

#endif // _ea452652_2d33_45ed_9feb_960f1d861041
