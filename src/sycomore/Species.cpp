#include "Species.h"

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

Species
::Species(Real R1, Real R2, Real D, Real R2_prime, Real delta_omega, Real w)
: R1(R1), R2(R2), D(D), R2_prime(R2_prime), delta_omega(delta_omega), w(w)
{
    // Nothing else.
}

Species
::Species(
    units::Frequency R1, units::Frequency R2,
    units::div<units::pow<units::Length, 2>, units::Time> D,
    units::Frequency R2_prime, Real delta_omega, Real w)
: R1(R1.convert_to(units::Hz)), R2(R2.convert_to(units::Hz)),
    D(D.convert_to(units::m*units::m/units::s)),
    R2_prime(R2_prime.convert_to(units::Hz)), delta_omega(delta_omega), w(w)
{
    // Nothing else.
}

}
