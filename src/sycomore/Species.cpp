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
    Quantity const & R1, Quantity const & R2, Quantity const & D,
    Quantity const & R2_prime, Quantity const & delta_omega, Real w)
: D(D.convert_to(units::m*units::m/units::s)),
    delta_omega(delta_omega.convert_to(units::rad/units::s)), w(w)
{
    this->R1 = (R1.dimensions==Time ? 1/R1 : R1).convert_to(units::Hz);
    this->R2 = (R2.dimensions==Time ? 1/R2 : R2).convert_to(units::Hz);
    this->R2_prime =
        (R2_prime.dimensions==Time ? 1/R2_prime : R2_prime).convert_to(
            units::Hz);
}

}
