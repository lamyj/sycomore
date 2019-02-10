#include "magnetization.h"

#include <cmath>
#include "sycomore/sycomore.h"

namespace sycomore
{

ComplexMagnetization const
ComplexMagnetization
::zero{{0,0}, 0, {0,0}};

ComplexMagnetization
::ComplexMagnetization()
: ComplexMagnetization(0, 0, 0)
{
    // Nothing else.
}

ComplexMagnetization
::ComplexMagnetization(Complex const & p, Real const & z, Complex const & m)
: p(p), m(m), z(z)
{
    // Nothing else.
}

bool
ComplexMagnetization
::operator==(ComplexMagnetization const & other) const
{
    return (this->p == other.p && this->z == other.z && this->m == other.m);
}

bool
ComplexMagnetization
::operator!=(ComplexMagnetization const & other) const
{
    return !this->operator==(other);
}

Magnetization::value_type transversal(Magnetization const & m)
{
    return std::sqrt(std::pow(m[0], 2.) + std::pow(m[1], 2.));
}

ComplexMagnetization as_complex_magnetization(Magnetization const & m)
{
    return ComplexMagnetization(
        Complex(m[0], m[1])/std::sqrt(2.),
        m[2],
        Complex(m[0], -m[1])/std::sqrt(2.));
}

Magnetization as_real_magnetization(ComplexMagnetization const & m)
{
    return Magnetization{
        ((m.p+m.m) / std::sqrt(2.)).real(),
        ((m.p-m.m) / std::sqrt(2.)).imag(),
        m.z};
}

}
