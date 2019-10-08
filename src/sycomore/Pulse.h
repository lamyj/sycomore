#ifndef _b46b179a_64e7_41db_9643_82051c4aa85a
#define _b46b179a_64e7_41db_9643_82051c4aa85a

#include <functional>
#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

class SYCOMORE_API Pulse
{
public:
    using RotationMatrix = Grid<Complex>;

    Pulse(Quantity const & angle, Quantity const & phase);

    Quantity const & get_angle() const;
    void set_angle(Quantity const & q);

    Quantity const & get_phase() const;
    void set_phase(Quantity const & q);

    /// @brief Return the rotation matrix for complex magnetization.
    RotationMatrix rotation_matrix() const;
private:
    /// @brief Flip angle.
    Quantity _angle;

    /// @brief Phase.
    Quantity _phase;
};

}

#endif // _b46b179a_64e7_41db_9643_82051c4aa85a
