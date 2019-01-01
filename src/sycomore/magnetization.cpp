#include "magnetization.h"

#include "sycomore/sycomore.h"

namespace sycomore
{

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
::ComplexMagnetization(
    Complex const & plus, Real const & zero, Complex const & minus)
: plus(plus), minus(minus), zero(zero)
{
    // Nothing else.
}

bool
ComplexMagnetization
::operator==(ComplexMagnetization const & other) const
{
    return (
        this->plus == other.plus
        && this->zero == other.zero
        && this->minus == other.minus);
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
        ((m.plus+m.minus) / std::sqrt(2.)).real(),
        ((m.plus-m.minus) / std::sqrt(2.)).imag(),
        m.zero};
}

}
