#ifndef _b46b179a_64e7_41db_9643_82051c4aa85a
#define _b46b179a_64e7_41db_9643_82051c4aa85a

#include <eigen3/Eigen/Dense>
#include "sycomore/sycomore.h"

namespace sycomore
{

class Pulse
{
public:
    using RotationMatrix = Eigen::Matrix<Complex, 3, 3>;

    /// @brief Flip angle in radians.
    Real angle;

    /// @brief Phase in radians.
    Real phase;

    /// @brief Return the rotation matrix for complex magnetization.
    RotationMatrix rotation_matrix() const;
};

}

#endif // _b46b179a_64e7_41db_9643_82051c4aa85a
