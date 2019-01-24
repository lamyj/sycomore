#include "Pulse.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/sycomore.h"

namespace sycomore
{

Pulse
::Pulse(double angle, double phase)
: angle(angle), phase(phase)
{
    // Nothing else.
}

Pulse
::Pulse(units::Angle angle, units::Angle phase)
: Pulse(angle.convert_to(units::rad), phase.convert_to(units::rad))
{
    // Nothing else.
}

Pulse::RotationMatrix
Pulse
::rotation_matrix() const
{
    Real const c_alpha = std::cos(this->angle);
    Real const s_alpha = std::sin(this->angle);
    Complex const e_i_phi{std::cos(this->phase), std::sin(this->phase)};

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

std::vector<sycomore::Pulse>
hard_pulse_approximation(
    Pulse const & pulse,
    std::function<Real(Real)> const & envelope,
    std::vector<Real> const & support)
{
    std::vector<Real> angles(support.size());
    std::transform(support.begin(), support.end(), angles.begin(), envelope);
    auto const sum = std::accumulate(angles.begin(), angles.end(), 0.);
    std::transform(
        angles.begin(), angles.end(), angles.begin(),
        [&](Real x) { return x*pulse.angle / sum; });

    std::vector<Pulse> pulses;
    for(auto && angle: angles)
    {
        pulses.emplace_back(angle, pulse.phase);
    }

    return pulses;
}

}
