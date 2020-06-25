#ifndef _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
#define _0bc5dc9b_ebb8_4139_bd22_f07f58e07314

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

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
        Quantity const & R2_prime=0*units::Hz,
        Quantity const & delta_omega=0*units::Hz, Real w=1);

    Species(
        Quantity const & R1, Quantity const & R2,
        Array<Quantity> const & D,
        Quantity const & R2_prime=0*units::Hz,
        Quantity const & delta_omega=0*units::Hz, Real w=1);

    /// @brief Return the longitudinal relaxation rate.
    Quantity const & get_R1() const;
    /// @brief Set the longitudinal relaxation rate or time.
    void set_R1(Quantity const & q);
    /// @brief Return the longitudinal relaxation time.
    Quantity const & get_T1() const;

    /// @brief Return the transversal relaxation rate.
    Quantity const & get_R2() const;
    /// @brief Set the transversal relaxation rate or time.
    void set_R2(Quantity const & q);
    /// @brief Return the transversal relaxation time.
    Quantity const & get_T2() const;

    /// @brief Return the diffusion tensor.
    Array<Quantity> const & get_D() const;
    /// @brief Set the diffusion coefficient (i.e. diagonal diffusion tensor).
    void set_D(Quantity const & q);
    /// @brief Set the diffusion tensor.
    void set_D(Array<Quantity> const & q);

    /**
     * @brief Return the part of the apparent transversal relaxation rate R2* 
     * attributed to the magnetic field inhomogeneity.
     */
    Quantity const & get_R2_prime() const;
    /**
     * @brief Set the part of the apparent transversal relaxation rate or time  
     * attributed to the magnetic field inhomogeneity.
     */
    void set_R2_prime(Quantity const & q);
    /**
     * @brief Return the part of the apparent transversal relaxation time T2* 
     * attributed to the magnetic field inhomogeneity.
     */
    Quantity const & get_T2_prime() const;

    /// @brief Return the frequency offset.
    Quantity const & get_delta_omega() const;
    /// @brief Set the frequency offset.
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
