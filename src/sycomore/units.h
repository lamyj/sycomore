#ifndef _e25f009a_96c8_4c52_97b5_de94a0752e6c
#define _e25f009a_96c8_4c52_97b5_de94a0752e6c

#include <cmath>
#include <ostream>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

namespace units
{

/**
 * @addtogroup BasicUnits
 * @{
 */

#define SYCOMORE_DECLARE_UNIT(dimensions, name, factor) \
    SYCOMORE_API extern Quantity const name; \
    SYCOMORE_API Quantity operator "" _##name(unsigned long long v); \
    SYCOMORE_API Quantity operator "" _##name(long double v);

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
//SYCOMORE_DECLARE_UNIT(Time, min, 60)
SYCOMORE_DECLARE_UNIT(Time, h, 3600)

SYCOMORE_DECLARE_UNITS(ElectricCurrent, A)
SYCOMORE_DECLARE_UNITS(ThermodynamicTemperature, K)
SYCOMORE_DECLARE_UNITS(AmountOfSubstance, mol)
SYCOMORE_DECLARE_UNITS(LuminousIntensity, cd)

/// @}

/**
 * @addtogroup DerivedUnits
 * @{
 */

SYCOMORE_DECLARE_UNITS(Angle, rad);
SYCOMORE_DECLARE_UNIT(Angle, deg, M_PI/180.);

SYCOMORE_DECLARE_UNITS(SolidAngle, sr);

SYCOMORE_DECLARE_UNITS(Frequency, Hz);
SYCOMORE_DECLARE_UNITS(Force, N);
SYCOMORE_DECLARE_UNITS(Pressure, Pa);
SYCOMORE_DECLARE_UNITS(Energy, J);
SYCOMORE_DECLARE_UNITS(Power, W);
SYCOMORE_DECLARE_UNITS(ElectricCharge, C);
SYCOMORE_DECLARE_UNITS(Voltage, V);
SYCOMORE_DECLARE_UNITS(Capacitance, F);
SYCOMORE_DECLARE_UNITS(Resistance, Ohm);
SYCOMORE_DECLARE_UNITS(ElectricalConductance, S);
SYCOMORE_DECLARE_UNITS(MagneticFlux, Wb);

// WARNING: _MT is a macro on Windows
#ifdef _WIN32
#pragma push_macro("_MT")
#undef _MT
#endif

SYCOMORE_DECLARE_UNITS(MagneticFluxDensity, T);

#ifdef _WIN32
#pragma pop_macro("_MT")
#endif

SYCOMORE_DECLARE_UNIT(MagneticFluxDensity, G, 1e-4);
SYCOMORE_DECLARE_UNITS(Inductance, H);
SYCOMORE_DECLARE_UNITS(LuminousFlux, lm);
SYCOMORE_DECLARE_UNITS(Illuminance, lx);
SYCOMORE_DECLARE_UNITS(Radioactivity, Bq);
SYCOMORE_DECLARE_UNITS(AbsorbedDose, Gy);
SYCOMORE_DECLARE_UNITS(EquivalentDose, Sv);
SYCOMORE_DECLARE_UNITS(CatalyticActivity, kat);

/// @}

}

}

//#include "units.txx"

#endif // _e25f009a_96c8_4c52_97b5_de94a0752e6c
