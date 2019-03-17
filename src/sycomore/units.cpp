#include "units.h"

#include <cmath>
#include <ostream>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"

namespace sycomore
{

namespace units
{

#define SYCOMORE_DEFINE_UNIT(Type, name) \
    Quantity operator "" _##name(unsigned long long v) { return double(v)*name; } \
    Quantity operator "" _##name(long double v) { return double(v)*name; }

#define SYCOMORE_DEFINE_UNITS(Type, name) \
    SYCOMORE_DEFINE_UNIT(Type, Y##name) \
    SYCOMORE_DEFINE_UNIT(Type, Z##name) \
    SYCOMORE_DEFINE_UNIT(Type, E##name) \
    SYCOMORE_DEFINE_UNIT(Type, P##name) \
    SYCOMORE_DEFINE_UNIT(Type, T##name) \
    SYCOMORE_DEFINE_UNIT(Type, G##name) \
    SYCOMORE_DEFINE_UNIT(Type, M##name) \
    SYCOMORE_DEFINE_UNIT(Type, k##name) \
    SYCOMORE_DEFINE_UNIT(Type, h##name) \
    SYCOMORE_DEFINE_UNIT(Type, da##name) \
    SYCOMORE_DEFINE_UNIT(Type, name) \
    SYCOMORE_DEFINE_UNIT(Type, d##name) \
    SYCOMORE_DEFINE_UNIT(Type, c##name) \
    SYCOMORE_DEFINE_UNIT(Type, m##name) \
    SYCOMORE_DEFINE_UNIT(Type, u##name) \
    SYCOMORE_DEFINE_UNIT(Type, n##name) \
    SYCOMORE_DEFINE_UNIT(Type, p##name) \
    SYCOMORE_DEFINE_UNIT(Type, f##name) \
    SYCOMORE_DEFINE_UNIT(Type, a##name) \
    SYCOMORE_DEFINE_UNIT(Type, z##name) \
    SYCOMORE_DEFINE_UNIT(Type, y##name)

SYCOMORE_DEFINE_UNITS(Length, m)
// WARNING: the base unit is kg, but since we only define the cast operators,
// this simpler form works.
SYCOMORE_DEFINE_UNITS(Mass, g)
SYCOMORE_DEFINE_UNITS(Time, s)
//SYCOMORE_DEFINE_UNIT(Time, min)
SYCOMORE_DEFINE_UNIT(Time, h)
SYCOMORE_DEFINE_UNITS(ElectricCurrent, A)
SYCOMORE_DEFINE_UNITS(ThermodynamicTemperature, K)
SYCOMORE_DEFINE_UNITS(AmountOfSubstance, mol)
SYCOMORE_DEFINE_UNITS(LuminousIntensity, cd)

SYCOMORE_DEFINE_UNITS(Angle, rad)
SYCOMORE_DEFINE_UNIT(Angle, deg)
SYCOMORE_DEFINE_UNITS(SolidAngle, sr)

SYCOMORE_DEFINE_UNITS(Frequency, Hz)
SYCOMORE_DEFINE_UNITS(Force, N)
SYCOMORE_DEFINE_UNITS(Pressure, Pa)
SYCOMORE_DEFINE_UNITS(Energy, J)
SYCOMORE_DEFINE_UNITS(Power, W)
SYCOMORE_DEFINE_UNITS(ElectricCharge, C)
SYCOMORE_DEFINE_UNITS(Voltage, V)
SYCOMORE_DEFINE_UNITS(Capacitance, F)
SYCOMORE_DEFINE_UNITS(Resistance, Ohm)
SYCOMORE_DEFINE_UNITS(ElectricalConductance, S)
SYCOMORE_DEFINE_UNITS(MagneticFlux, Wb)
SYCOMORE_DEFINE_UNITS(MagneticFluxDensity, T)
SYCOMORE_DEFINE_UNIT(MagneticFluxDensity, G)
SYCOMORE_DEFINE_UNITS(Inductance, H)
SYCOMORE_DEFINE_UNITS(LuminousFlux, lm)
SYCOMORE_DEFINE_UNITS(Illuminance, lx)
SYCOMORE_DEFINE_UNITS(Radioactivity, Bq)
SYCOMORE_DEFINE_UNITS(AbsorbedDose, Gy)
SYCOMORE_DEFINE_UNITS(EquivalentDose, Sv)
SYCOMORE_DEFINE_UNITS(CatalyticActivity, kat)

}

}
