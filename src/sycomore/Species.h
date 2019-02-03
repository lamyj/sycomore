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

    /// @brief Relative frequency in rad/s.
    Real delta_omega;

    /// @brief Relative weight
    Real w;

    Species(
        Real R1, Real R2,
        Real D=0, Real R2_prime=0, Real delta_omega=0, Real w=1);

    Species(
        Quantity const & R1, Quantity const & R2,
        Quantity const & D={0, Diffusion}, Quantity const & R2_prime={0, Frequency},
        Quantity const & delta_omega={0, AngularFrequency}, Real w=1);
};

}

#endif // _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
