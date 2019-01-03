#ifndef _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
#define _0bc5dc9b_ebb8_4139_bd22_f07f58e07314

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

class Species
{
public:
    /// @brief R1 relaxation rate in s.
    Real R1;

    /// @brief R2 relaxation rate in s.
    Real R2;

    /// @brief Diffusivity in m^2 / s.
    Real D;

    /// @brief R2' relaxation rate in s.
    Real R2_prime;

    /// @brief Relative frequency
    Real delta_omega;

    /// @brief Relative weight
    Real w;

    Species(
        Real R1, Real R2,
        Real D=0, Real R2_prime=0, Real delta_omega=0, Real w=1);

    Species(
        units::Frequency R1, units::Frequency R2,
        units::div<units::pow<units::Length, 2>, units::Time> D=
            units::div<units::pow<units::Length, 2>, units::Time>(0.),
        units::Frequency R2_prime=units::Frequency(0.), Real delta_omega=0,
        Real w=1);

    Species(
        units::Time T1, units::Time T2,
        units::div<units::pow<units::Length, 2>, units::Time> D=
            units::div<units::pow<units::Length, 2>, units::Time>(0.),
        units::Time T2_prime=units::Time(INFINITY), Real delta_omega=0,
        Real w=1);
};

}

#endif // _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
