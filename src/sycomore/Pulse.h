#ifndef _b46b179a_64e7_41db_9643_82051c4aa85a
#define _b46b179a_64e7_41db_9643_82051c4aa85a

#include <functional>
#include <vector>

#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

namespace sycomore
{

class SYCOMORE_API Pulse
{
public:
    Pulse(Quantity const & angle, Quantity const & phase=0*units::rad);

    Quantity const & get_angle() const;
    void set_angle(Quantity const & q);

    Quantity const & get_phase() const;
    void set_phase(Quantity const & q);
private:
    /// @brief Flip angle.
    Quantity _angle;

    /// @brief Phase.
    Quantity _phase;
};

}

#endif // _b46b179a_64e7_41db_9643_82051c4aa85a
