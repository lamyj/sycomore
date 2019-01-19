#include "magnetization.h"

#include <cmath>
#include "sycomore/sycomore.h"

namespace sycomore
{

ComplexMagnetization const
ComplexMagnetization
::zero{{0,0}, 0, {0,0}};

Real
Magnetization
::transversal() const
{
    return std::sqrt(std::pow(this->x, Real(2.)) + std::pow(this->y, Real(2.)));
}

bool
Magnetization
::operator==(Magnetization const & other) const
{
    return this->x == other.x && this->y == other.y && this->z == other.z;
}

bool
Magnetization
::operator!=(Magnetization const & other) const
{
    return !this->operator==(other);
}

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

ComplexMagnetization as_complex_magnetization(Magnetization const & m)
{
    return ComplexMagnetization(
        Complex(m.x, m.y)/std::sqrt(2.),
        m.z,
        Complex(m.x, -m.y)/std::sqrt(2.));
}

Magnetization as_real_magnetization(ComplexMagnetization const & m)
{
    return Magnetization{
        ((m.p+m.m) / std::sqrt(2.)).real(),
        ((m.p-m.m) / std::sqrt(2.)).imag(),
        m.z};
}

}
