#ifndef _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
#define _0bc5dc9b_ebb8_4139_bd22_f07f58e07314

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

class SYCOMORE_API Species
{
public:
    /// @brief Relative weight
    Real w;

    Species(Quantity const & R1, Quantity const & R2);

    Species(
        Quantity const & R1, Quantity const & R2,
        Quantity const & D,
        Quantity const & R2_prime={0, Frequency},
        Quantity const & delta_omega={0, AngularFrequency}, Real w=1);

    Species(
        Quantity const & R1, Quantity const & R2,
        Array<Quantity> const & D,
        Quantity const & R2_prime={0, Frequency},
        Quantity const & delta_omega={0, AngularFrequency}, Real w=1);

    Quantity const & get_R1() const;
    void set_R1(Quantity const & q);
    Quantity const & get_T1() const;

    Quantity const & get_R2() const;
    void set_R2(Quantity const & q);
    Quantity const & get_T2() const;

    Array<Quantity> const & get_D() const;
    void set_D(Quantity const & q);
    void set_D(Array<Quantity> const & q);

    Quantity const & get_R2_prime() const;
    void set_R2_prime(Quantity const & q);
    Quantity const & get_T2_prime() const;

    Quantity const & get_delta_omega() const;
    void set_delta_omega(Quantity const & q);

private:
    /// @brief R1 relaxation rate.
    Quantity _R1;

    /// @brief R1 relaxation time.
    Quantity _T1;

    /// @brief R2 relaxation rate.
    Quantity _R2;

    /// @brief R2 relaxation time.
    Quantity _T2;

    /// @brief Diffusivity.
    Array<Quantity> _D;

    /// @brief R2' relaxation rate.
    Quantity _R2_prime;

    /// @brief R2' relaxation time.
    Quantity _T2_prime;

    /// @brief Relative frequency.
    Quantity _delta_omega;
};

}

#endif // _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
