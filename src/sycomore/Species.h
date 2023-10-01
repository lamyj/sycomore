#ifndef _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
#define _0bc5dc9b_ebb8_4139_bd22_f07f58e07314

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

/// @brief Species described by its NMR parameters
class Species
{
public:
    /// @brief Create a species given relaxation rates or times
    Species(Quantity const & R1, Quantity const & R2);
    
    /**
     * @brief Create a species given relaxation rates or times, isotropic
     * diffusion coefficient, off-resonance and weight.
     */
    Species(
        Quantity const & R1, Quantity const & R2, Quantity const & D,
        Quantity const & delta_omega=0*units::Hz);

    /**
     * @brief Create a species given relaxation rates or times, diffusion
     * tensor, off-resonance and weight.
     */
    Species(
        Quantity const & R1, Quantity const & R2, Matrix3x3Q const & D,
        Quantity const & delta_omega=0*units::Hz);

    /// @brief Return the longitudinal relaxation rate.
    Quantity const & R1() const;
    /// @brief Set the longitudinal relaxation rate or time.
    void set_R1(Quantity const & q);
    /// @brief Return the longitudinal relaxation time.
    Quantity const & T1() const;
    /// @brief Set the longitudinal relaxation rate or time.
    void set_T1(Quantity const & q);
    
    /// @brief Return the transversal relaxation rate.
    Quantity const & R2() const;
    /// @brief Set the transversal relaxation rate or time.
    void set_R2(Quantity const & q);
    /// @brief Return the transversal relaxation time.
    Quantity const & T2() const;
    /// @brief Set the transversal relaxation rate or time.
    void set_T2(Quantity const & q);

    /// @brief Return the diffusion tensor.
    Matrix3x3Q const & D() const;
    /// @brief Set the diffusion coefficient (i.e. diagonal diffusion tensor).
    void set_D(Quantity const & q);
    /// @brief Set the diffusion tensor.
    void set_D(Matrix3x3Q const & q);

    /// @brief Return the frequency offset.
    Quantity const & delta_omega() const;
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
    Matrix3x3Q _D;

    /// @brief Relative frequency.
    Quantity _delta_omega;
};

}

#endif // _0bc5dc9b_ebb8_4139_bd22_f07f58e07314
