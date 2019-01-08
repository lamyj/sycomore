#ifndef _84b1fab9_d85f_41d4_9251_3e1b61073ec5
#define _84b1fab9_d85f_41d4_9251_3e1b61073ec5

#include "sycomore/sycomore.h"

namespace sycomore
{

struct Magnetization
{
    Real x, y, z;

    /// @brief Return the transversal magnetization.
    Real transversal() const;

    bool operator==(Magnetization const & other) const;
    bool operator!=(Magnetization const & other) const;
};

struct ComplexMagnetization
{
    static ComplexMagnetization const zero;

    Complex p;
    Real z;
    Complex m;

    ComplexMagnetization();

    ComplexMagnetization(Complex const & p, Real const & z, Complex const & m);

    bool operator==(ComplexMagnetization const & other) const;
    bool operator!=(ComplexMagnetization const & other) const;
};

ComplexMagnetization as_complex_magnetization(Magnetization const & m);
Magnetization as_real_magnetization(ComplexMagnetization const & m);

}

#endif // _84b1fab9_d85f_41d4_9251_3e1b61073ec5
