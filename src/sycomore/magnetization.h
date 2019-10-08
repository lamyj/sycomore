#ifndef _84b1fab9_d85f_41d4_9251_3e1b61073ec5
#define _84b1fab9_d85f_41d4_9251_3e1b61073ec5

#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

struct SYCOMORE_API ComplexMagnetization
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
using Magnetization = Array<Real>;

SYCOMORE_API Magnetization::value_type transversal(Magnetization const & m);
SYCOMORE_API ComplexMagnetization as_complex_magnetization(Magnetization const & m);
SYCOMORE_API Magnetization as_real_magnetization(ComplexMagnetization const & m);

}

#endif // _84b1fab9_d85f_41d4_9251_3e1b61073ec5
