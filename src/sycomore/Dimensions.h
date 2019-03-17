#ifndef _bad9a44e_bb0e_403b_be26_3ca5d92c4d4a
#define _bad9a44e_bb0e_403b_be26_3ca5d92c4d4a

#include <ostream>

namespace sycomore
{

class Dimensions
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

Dimensions operator*(Dimensions l, Dimensions const & r);
Dimensions operator/(Dimensions l, Dimensions const & r);

std::ostream & operator<<(std::ostream & stream, Dimensions const & d);

}

namespace std
{

sycomore::Dimensions pow(sycomore::Dimensions const & d, double s);

}

namespace sycomore
{

/**
 * @addtogroup BasicDimensions.
 * @{
 */
Dimensions const Length{1, 0, 0, 0, 0, 0, 0};
Dimensions const Mass{0, 1, 0, 0, 0, 0, 0};
Dimensions const Time{0, 0, 1, 0, 0, 0, 0};
Dimensions const ElectricCurrent{0, 0, 0, 1, 0, 0, 0};
Dimensions const ThermodynamicTemperature{0, 0, 0, 0, 1, 0, 0};
Dimensions const AmountOfSubstance{0, 0, 0, 0, 0, 1, 0};
Dimensions const LuminousIntensity{0, 0, 0, 0, 0, 0, 1};
/// @}

/**
 * @addtogroup DerivedDimensions
 * @{
 */
Dimensions const Surface = std::pow(Length, 2);
Dimensions const Volume = std::pow(Length, 3);

Dimensions const Velocity = Length/Time;
Dimensions const Acceleration = Velocity/Time;

Dimensions const Angle = Length/Length;
Dimensions const SolidAngle = Surface/Surface;

Dimensions const Frequency = std::pow(Time, -1);
Dimensions const Force = Mass*Length*std::pow(Time, -2);
Dimensions const Pressure = Force/Surface;
Dimensions const Energy = Force*Length;
Dimensions const Power = Energy/Time;
Dimensions const ElectricCharge = Time*ElectricCurrent;
Dimensions const Voltage = Power/ElectricCurrent;
Dimensions const Capacitance = ElectricCharge/Voltage;
Dimensions const Resistance = Voltage/ElectricCurrent;
Dimensions const ElectricalConductance = std::pow(Resistance, -1);
Dimensions const MagneticFlux = Energy/ElectricCurrent;
Dimensions const MagneticFluxDensity = MagneticFlux/Surface;
Dimensions const Inductance = MagneticFlux/ElectricCurrent;
Dimensions const LuminousFlux = LuminousIntensity/SolidAngle;
Dimensions const Illuminance = LuminousFlux/Surface;
Dimensions const Radioactivity = std::pow(Time, -1);
Dimensions const AbsorbedDose = Energy/Mass;
Dimensions const EquivalentDose = Energy/Mass;
Dimensions const CatalyticActivity = AmountOfSubstance/Time;

Dimensions const AngularFrequency = Angle/Time;

/// @}

}

#endif // _bad9a44e_bb0e_403b_be26_3ca5d92c4d4a
