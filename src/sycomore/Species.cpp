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
    units::Frequency R1, units::Frequency R2, Diffusion D,
    units::Frequency R2_prime, units::AngularFrequency delta_omega, Real w)
: R1(R1.convert_to(units::Hz)), R2(R2.convert_to(units::Hz)),
    D(D.convert_to(units::m*units::m/units::s)),
    R2_prime(R2_prime.convert_to(units::Hz)),
    delta_omega(delta_omega.convert_to(units::rad/units::s)), w(w)
{
    // Nothing else.
}

Species
::Species(
    units::Time T1, units::Time T2, Diffusion D,
    units::Time T2_prime, units::AngularFrequency delta_omega, Real w)
: R1((1./T1).convert_to(units::Hz)), R2((1./T2).convert_to(units::Hz)),
    D(D.convert_to(units::m*units::m/units::s)),
    R2_prime((1./T2_prime).convert_to(units::Hz)),
    delta_omega(delta_omega.convert_to(units::rad/units::s)), w(w)
{
    // Nothing else.
}

}
