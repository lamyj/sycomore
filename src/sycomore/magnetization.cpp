#include "sycomore/magnetization.h"

#include <eigen3/Eigen/Dense>
#include "sycomore/sycomore.h"

namespace sycomore
{

ComplexMagnetization as_complex_magnetization(Magnetization const & m)
{
    return ComplexMagnetization(
        Complex(m(0), m(1))/std::sqrt(2.),
        m(2),
        Complex(m(0), -m(1))/std::sqrt(2.));
}

Magnetization as_real_magnetization(ComplexMagnetization const & m)
{
    return Magnetization(
        ((m(0)+m(2)) / std::sqrt(2.)).real(),
        ((m(0)-m(2)) / std::sqrt(2.)).imag(),
        m(1).real());
}

}
