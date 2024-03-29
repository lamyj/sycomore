#ifndef _ea452652_2d33_45ed_9feb_960f1d861041
#define _ea452652_2d33_45ed_9feb_960f1d861041

#include <functional>
#include <string>
#include <vector>

#include "sycomore/Pulse.h"
#include "sycomore/Quantity.h"
#include "sycomore/units.h"

namespace sycomore
{

/// @brief Small tip angle approximation of a shaped pulse
class HardPulseApproximation
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
    std::vector<Pulse> const & pulses() const;
    
    /// @brief Return the duration of a hard pulse
    Quantity const & duration() const;
    
    /// @brief Return the phase of the shaped pulse
    Quantity const & phase() const;
    
    /// @brief Set the phase of the shaped pulse
    void set_phase(Quantity const & phase);

private:
    std::vector<Pulse> _pulses;
    Quantity _duration;
};

/// @addtogroup HardPulseApproximationEnvelopes
/// @{

/// @brief Create an apodized sinc envelope
HardPulseApproximation::Envelope apodized_sinc_envelope(
    Quantity const & t0, unsigned int N, Real alpha);
/// @brief Create an Hann-apodized sinc envelope
HardPulseApproximation::Envelope hann_sinc_envelope(
    Quantity const & t0, unsigned int N);
/// @brief Create an Hamming-apodized sinc envelope
HardPulseApproximation::Envelope hamming_sinc_envelope(
    Quantity const & t0, unsigned int N);
/// @brief Create a sinc envelope
HardPulseApproximation::Envelope sinc_envelope(Quantity const & t0);

/// @}

}

#endif // _ea452652_2d33_45ed_9feb_960f1d861041
