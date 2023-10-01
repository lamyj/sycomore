#ifndef _b46b179a_64e7_41db_9643_82051c4aa85a
#define _b46b179a_64e7_41db_9643_82051c4aa85a

#include <functional>
#include <vector>

#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

/// @brief RF pulse
class Pulse
{
public:
    /// @brief Creator
    Pulse(Quantity const & angle, Quantity const & phase=0*units::rad);
    
    /// @brief Return the flip angle of the pulse
    Quantity const & angle() const;
    /// @brief Set the flip angle of the pulse
    void set_angle(Quantity const & q);
    
    /// @brief Return the phase of the pulse
    Quantity const & phase() const;
    /// @brief Set the phase of the pulse
    void set_phase(Quantity const & q);
private:
    /// @brief Flip angle.
    Quantity _angle;

    /// @brief Phase.
    Quantity _phase;
};

}

#endif // _b46b179a_64e7_41db_9643_82051c4aa85a
