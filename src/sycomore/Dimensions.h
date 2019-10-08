#ifndef _bad9a44e_bb0e_403b_be26_3ca5d92c4d4a
#define _bad9a44e_bb0e_403b_be26_3ca5d92c4d4a

#include <ostream>

#include "sycomore/sycomore_api.h"

namespace sycomore
{

class SYCOMORE_API Dimensions
{
public:
    double length;
    double mass;
    double time;
    double electric_current;
    double thermodynamic_temperature;
    double amount_of_substance;
    double luminous_intensity;

    Dimensions(
        double length=0, double mass=0, double time=0, double electric_current=0,
        double thermodynamic_temperature=0, double amount_of_substance=0,
        double luminous_intensity=0);

    bool operator==(Dimensions const & other) const;
    bool operator!=(Dimensions const & other) const;

    Dimensions & operator*=(Dimensions const & other);
    Dimensions & operator/=(Dimensions const & other);
};

SYCOMORE_API extern Dimensions operator*(Dimensions l, Dimensions const & r);
SYCOMORE_API extern Dimensions operator/(Dimensions l, Dimensions const & r);

SYCOMORE_API std::ostream & operator<<(std::ostream & stream, Dimensions const & d);

}

namespace std
{

SYCOMORE_API sycomore::Dimensions pow(sycomore::Dimensions const & d, double s);

}

namespace sycomore
{

/**
 * @addtogroup BasicDimensions.
 * @{
 */
SYCOMORE_API extern Dimensions const Length;
SYCOMORE_API extern Dimensions const Mass;
SYCOMORE_API extern Dimensions const Time;
SYCOMORE_API extern Dimensions const ElectricCurrent;
SYCOMORE_API extern Dimensions const ThermodynamicTemperature;
SYCOMORE_API extern Dimensions const AmountOfSubstance;
SYCOMORE_API extern Dimensions const LuminousIntensity;
/// @}

/**
 * @addtogroup DerivedDimensions
 * @{
 */
SYCOMORE_API extern Dimensions const Surface;
SYCOMORE_API extern Dimensions const Volume;

SYCOMORE_API extern Dimensions const Velocity;
SYCOMORE_API extern Dimensions const Acceleration;

SYCOMORE_API extern Dimensions const Angle;
SYCOMORE_API extern Dimensions const SolidAngle;

SYCOMORE_API extern Dimensions const Frequency;
SYCOMORE_API extern Dimensions const Force;
SYCOMORE_API extern Dimensions const Pressure;
SYCOMORE_API extern Dimensions const Energy;
SYCOMORE_API extern Dimensions const Power;
SYCOMORE_API extern Dimensions const ElectricCharge;
SYCOMORE_API extern Dimensions const Voltage;
SYCOMORE_API extern Dimensions const Capacitance;
SYCOMORE_API extern Dimensions const Resistance;
SYCOMORE_API extern Dimensions const ElectricalConductance;
SYCOMORE_API extern Dimensions const MagneticFlux;
SYCOMORE_API extern Dimensions const MagneticFluxDensity;
SYCOMORE_API extern Dimensions const Inductance;
SYCOMORE_API extern Dimensions const LuminousFlux;
SYCOMORE_API extern Dimensions const Illuminance;
SYCOMORE_API extern Dimensions const Radioactivity;
SYCOMORE_API extern Dimensions const AbsorbedDose;
SYCOMORE_API extern Dimensions const EquivalentDose;
SYCOMORE_API extern Dimensions const CatalyticActivity;

SYCOMORE_API extern Dimensions const AngularFrequency;

/// @}

}

#endif // _bad9a44e_bb0e_403b_be26_3ca5d92c4d4a
