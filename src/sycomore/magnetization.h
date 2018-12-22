#ifndef _84b1fab9_d85f_41d4_9251_3e1b61073ec5
#define _84b1fab9_d85f_41d4_9251_3e1b61073ec5

#include <eigen3/Eigen/Dense>
#include "sycomore/sycomore.h"

namespace sycomore
{

using Real = double;
using Complex = std::complex<Real>;

using Magnetization = Eigen::Vector3d;
using ComplexMagnetization = Eigen::Matrix<Complex, 3, 1>;

ComplexMagnetization as_complex_magnetization(Magnetization const & m);
Magnetization as_real_magnetization(ComplexMagnetization const & m);

}

#endif // _84b1fab9_d85f_41d4_9251_3e1b61073ec5
