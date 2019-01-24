#ifndef _b46b179a_64e7_41db_9643_82051c4aa85a
#define _b46b179a_64e7_41db_9643_82051c4aa85a

#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

class Pulse
{
public:
    using RotationMatrix = Grid<Complex>;

    /// @brief Flip angle in radians.
    Real angle;

    /// @brief Phase in radians.
    Real phase;

    Pulse(double angle, double phase);

    Pulse(units::Angle angle, units::Angle phase);

    /// @brief Return the rotation matrix for complex magnetization.
    RotationMatrix rotation_matrix() const;
};

std::vector<Pulse>
hard_pulse_approximation(
    Pulse const & pulse, std::function<Real(Real)> const & envelope,
    std::vector<Real> const & support);

}

#endif // _b46b179a_64e7_41db_9643_82051c4aa85a
