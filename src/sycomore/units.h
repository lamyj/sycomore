#ifndef _e25f009a_96c8_4c52_97b5_de94a0752e6c
#define _e25f009a_96c8_4c52_97b5_de94a0752e6c

#include <cmath>
#include <ostream>

namespace sycomore
{

namespace units
{

template<int L, int M, int T, int I, int Theta, int N, int J>
class Unit
{
public:
    using Self = Unit<L,M,T,I,Theta,N,J>;

    double value;

    Unit();
    explicit Unit(double value);

    Unit(Self const &) = default;
    Unit(Self &&) = default;
    Self & operator=(Self const &) = default;
    Self & operator=(Self &&) = default;
    ~Unit() = default;

    double convert_to(Unit dest) const;
};

/******************************************************************************
 ******************************* Combining units ******************************
 ******************************************************************************/

template<typename, typename> struct _mult;
template<
    int L1, int M1, int T1, int I1, int Theta1, int N1, int J1,
    int L2, int M2, int T2, int I2, int Theta2, int N2, int J2>
struct _mult<Unit<L1,M1,T1,I1,Theta1,N1,J1>, Unit<L2,M2,T2,I2,Theta2,N2,J2>>
{
    using type = Unit<L1+L2,M1+M2,T1+T2,I1+I2,Theta1+Theta2,N1+N2,J1+J2>;
};
template<typename U1, typename U2> using mult = typename _mult<U1, U2>::type;

template<typename, typename> struct _div;
template<
    int L1, int M1, int T1, int I1, int Theta1, int N1, int J1,
    int L2, int M2, int T2, int I2, int Theta2, int N2, int J2>
struct _div<Unit<L1,M1,T1,I1,Theta1,N1,J1>, Unit<L2,M2,T2,I2,Theta2,N2,J2>>
{
    using type = Unit<L1-L2,M1-M2,T1-T2,I1-I2,Theta1-Theta2,N1-N2,J1-J2>;
};
template<typename U1, typename U2> using div = typename _div<U1, U2>::type;

template<typename, int> struct _pow;
template<
    int L1, int M1, int T1, int I1, int Theta1, int N1, int J1,
    int Power>
struct _pow<Unit<L1,M1,T1,I1,Theta1,N1,J1>, Power>
{
    using type = Unit<L1*Power,M1*Power,T1*Power,I1*Power,Theta1*Power,N1*Power,J1*Power>;
};
template<typename U, int P> using pow = typename _pow<U, P>::type;

/******************************************************************************
 ******************************* Unit operators *******************************
 ******************************************************************************/

template<int ... Args>
Unit<Args ...>
operator-(Unit<Args ...> const & x);

template<int ... Args>
Unit<Args ...>
operator+(Unit<Args ...> const & x);

template<int ... Args>
Unit<Args ...>
operator+(Unit<Args ...> const & x, Unit<Args ...> const & y);

template<int ... Args>
Unit<Args ...>
operator-(Unit<Args ...> const & x, Unit<Args ...> const & y);

template<int ... Args>
Unit<Args ...>
operator*(double const & scalar, Unit<Args ...> const & unit);

template<int ... Args>
Unit<Args ...>
operator*(Unit<Args ...> const & unit, double const & scalar);

template<int ... Args>
pow<Unit<Args ...>, -1>
operator/(double const & scalar, Unit<Args ...> const & unit);

template<int ... Args>
Unit<Args ...>
operator/(Unit<Args ...> const & unit, double const & scalar);

template<typename Unit1, typename Unit2>
mult<Unit1, Unit2>
operator*(Unit1 const & unit_1, Unit2 const & unit_2);

template<typename Unit1, typename Unit2>
div<Unit1, Unit2>
operator/(Unit1 const & unit_1, Unit2 const & unit_2);

template<int L, int M, int T, int I, int Theta, int N, int J>
std::ostream &
operator<<(std::ostream & stream, Unit<L,M,T,I,Theta,N,J> const & u);

/******************************************************************************
 ********************************* Basic units ********************************
 ******************************************************************************/

using Length = Unit<1, 0, 0, 0, 0, 0, 0>;
using Mass = Unit<0, 1, 0, 0, 0, 0, 0>;
using Time = Unit<0, 0, 1, 0, 0, 0, 0>;
using ElectricCurrent = Unit<0, 0, 0, 1, 0, 0, 0>;
using ThermodynamicTemperature = Unit<0, 0, 0, 0, 1, 0, 0>;
using AmountOfSubstance = Unit<0, 0, 0, 0, 0, 1, 0>;
using LuminousIntensity = Unit<0, 0, 0, 0, 0, 0, 1>;

#define SYCOMORE_DECLARE_UNIT(Type, name, factor) \
    Type const name(factor); \
    Type operator "" _##name(unsigned long long v); \
    Type operator "" _##name(long double v);

#define SYCOMORE_DECLARE_UNITS(Type, name) \
    SYCOMORE_DECLARE_UNIT(Type, Y##name, 1e24) \
    SYCOMORE_DECLARE_UNIT(Type, Z##name, 1e21) \
    SYCOMORE_DECLARE_UNIT(Type, E##name, 1e18) \
    SYCOMORE_DECLARE_UNIT(Type, P##name, 1e15) \
    SYCOMORE_DECLARE_UNIT(Type, T##name, 1e12) \
    SYCOMORE_DECLARE_UNIT(Type, G##name, 1e9) \
    SYCOMORE_DECLARE_UNIT(Type, M##name, 1e6) \
    SYCOMORE_DECLARE_UNIT(Type, k##name, 1e3) \
    SYCOMORE_DECLARE_UNIT(Type, h##name, 1e2) \
    SYCOMORE_DECLARE_UNIT(Type, da##name, 1e1) \
    SYCOMORE_DECLARE_UNIT(Type, name, 1) \
    SYCOMORE_DECLARE_UNIT(Type, d##name, 1e-1) \
    SYCOMORE_DECLARE_UNIT(Type, c##name, 1e-2) \
    SYCOMORE_DECLARE_UNIT(Type, m##name, 1e-3) \
    SYCOMORE_DECLARE_UNIT(Type, u##name, 1e-6) \
    SYCOMORE_DECLARE_UNIT(Type, n##name, 1e-9) \
    SYCOMORE_DECLARE_UNIT(Type, p##name, 1e-12) \
    SYCOMORE_DECLARE_UNIT(Type, f##name, 1e-15) \
    SYCOMORE_DECLARE_UNIT(Type, a##name, 1e-18) \
    SYCOMORE_DECLARE_UNIT(Type, z##name, 1e-21) \
    SYCOMORE_DECLARE_UNIT(Type, y##name, 1e-24)

SYCOMORE_DECLARE_UNITS(Length, m)

SYCOMORE_DECLARE_UNIT(Mass, Yg, 1e21)
SYCOMORE_DECLARE_UNIT(Mass, Zg, 1e18)
SYCOMORE_DECLARE_UNIT(Mass, Eg, 1e15)
SYCOMORE_DECLARE_UNIT(Mass, Pg, 1e12)
SYCOMORE_DECLARE_UNIT(Mass, Tg, 1e9)
SYCOMORE_DECLARE_UNIT(Mass, Gg, 1e6)
SYCOMORE_DECLARE_UNIT(Mass, Mg, 1e3)
SYCOMORE_DECLARE_UNIT(Mass, kg, 1)
SYCOMORE_DECLARE_UNIT(Mass, hg, 1e-1)
SYCOMORE_DECLARE_UNIT(Mass, dag, 1e-2)
SYCOMORE_DECLARE_UNIT(Mass, g, 1e-3)
SYCOMORE_DECLARE_UNIT(Mass, dg, 1e-4)
SYCOMORE_DECLARE_UNIT(Mass, cg, 1e-5)
SYCOMORE_DECLARE_UNIT(Mass, mg, 1e-6)
SYCOMORE_DECLARE_UNIT(Mass, ug, 1e-9)
SYCOMORE_DECLARE_UNIT(Mass, ng, 1e-12)
SYCOMORE_DECLARE_UNIT(Mass, pg, 1e-15)
SYCOMORE_DECLARE_UNIT(Mass, fg, 1e-18)
SYCOMORE_DECLARE_UNIT(Mass, ag, 1e-21)
SYCOMORE_DECLARE_UNIT(Mass, zg, 1e-24)
SYCOMORE_DECLARE_UNIT(Mass, yg, 1e-27)

SYCOMORE_DECLARE_UNITS(Time, s)
SYCOMORE_DECLARE_UNIT(Time, min, 60)
SYCOMORE_DECLARE_UNIT(Time, h, 3600)

SYCOMORE_DECLARE_UNITS(ElectricCurrent, A)
SYCOMORE_DECLARE_UNITS(ThermodynamicTemperature, K)
SYCOMORE_DECLARE_UNITS(AmountOfSubstance, mol)
SYCOMORE_DECLARE_UNITS(LuminousIntensity, cd)

/******************************************************************************
 ******************************** Derived units *******************************
 ******************************************************************************/

using Surface = pow<Length, 2>;
using Volume = pow<Length, 3>;

using Velocity = div<Length, Time>;
using Acceleration = div<Length, pow<Time, 2>>;

#define SYCOMORE_DECLARE_DERIVED_UNIT(Type, name, ...) \
    using Type = __VA_ARGS__; \
    SYCOMORE_DECLARE_UNITS(Type, name);

SYCOMORE_DECLARE_DERIVED_UNIT(Angle, rad, div<Length, Length>);
// WARNING M_PI is not C++-standard, but is POSIX-standard.
// On Windows, define _USE_MATH_DEFINES
SYCOMORE_DECLARE_UNIT(Angle, deg, M_PI/180.);

SYCOMORE_DECLARE_DERIVED_UNIT(SolidAngle, sr, div<Surface, Surface>);

SYCOMORE_DECLARE_DERIVED_UNIT(Frequency, Hz, pow<Time, -1>);
SYCOMORE_DECLARE_DERIVED_UNIT(Force, N, div<mult<Mass, Length>, pow<Time, 2>>);
SYCOMORE_DECLARE_DERIVED_UNIT(Pressure, Pa, div<Force, Surface>);
SYCOMORE_DECLARE_DERIVED_UNIT(Energy, J, mult<Force, Length>);
SYCOMORE_DECLARE_DERIVED_UNIT(Power, W, div<Energy, Time>);
SYCOMORE_DECLARE_DERIVED_UNIT(ElectricCharge, C, mult<Time, ElectricCurrent>);
SYCOMORE_DECLARE_DERIVED_UNIT(Voltage, V, mult<Power, ElectricCurrent>);
SYCOMORE_DECLARE_DERIVED_UNIT(Capacitance, F, div<ElectricCharge, Voltage>);
SYCOMORE_DECLARE_DERIVED_UNIT(Resistance, Ohm, div<Voltage, ElectricCurrent>);
SYCOMORE_DECLARE_DERIVED_UNIT(ElectricalConductance, S, pow<Resistance, -1>);
SYCOMORE_DECLARE_DERIVED_UNIT(MagneticFlux, Wb, mult<Voltage, Time>);
SYCOMORE_DECLARE_DERIVED_UNIT(MagneticFluxDensity, T, div<MagneticFlux, Surface>);
SYCOMORE_DECLARE_DERIVED_UNIT(Inductance, H, div<MagneticFlux, ElectricCurrent>);
SYCOMORE_DECLARE_DERIVED_UNIT(LuminousFlux, lm, mult<LuminousIntensity, SolidAngle>);
SYCOMORE_DECLARE_DERIVED_UNIT(Illuminance, lx, mult<Surface, LuminousFlux>);
SYCOMORE_DECLARE_DERIVED_UNIT(Radioactivity, Bq, pow<Time, -1>);
SYCOMORE_DECLARE_DERIVED_UNIT(AbsorbedDose, Gy, div<Energy, Mass>);
SYCOMORE_DECLARE_DERIVED_UNIT(EquivalentDose, Sv, div<Energy, Mass>);
SYCOMORE_DECLARE_DERIVED_UNIT(CatalyticActivity, kat, div<AmountOfSubstance, Time>);

}

}

#include "units.txx"

#endif // _e25f009a_96c8_4c52_97b5_de94a0752e6c
