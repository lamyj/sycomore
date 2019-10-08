#include "Dimensions.h"

#include <ostream>
#include <string>

namespace sycomore
{

Dimensions
::Dimensions(
    double length, double mass, double time, double electric_current,
    double thermodynamic_temperature, double amount_of_substance,
    double luminous_intensity)
: length(length), mass(mass), time(time), electric_current(electric_current),
    thermodynamic_temperature(thermodynamic_temperature),
    amount_of_substance(amount_of_substance),
    luminous_intensity(luminous_intensity)
{
    // Nothing else.
}

bool
Dimensions
::operator==(Dimensions const & other) const
{
    return (
        this->length == other.length
        && this->mass == other.mass
        && this->time == other.time
        && this->electric_current == other.electric_current
        && this->thermodynamic_temperature == other.thermodynamic_temperature
        && this->amount_of_substance == other.amount_of_substance
        && this->luminous_intensity == other.luminous_intensity);
}

Dimensions &
Dimensions
::operator*=(Dimensions const & other)
{
    this->length += other.length;
    this->mass += other.mass;
    this->time += other.time;
    this->electric_current += other.electric_current;
    this->thermodynamic_temperature += other.thermodynamic_temperature;
    this->amount_of_substance += other.amount_of_substance;
    this->luminous_intensity += other.luminous_intensity;

    return *this;
}

Dimensions &
Dimensions
::operator/=(Dimensions const & other)
{
    this->length -= other.length;
    this->mass -= other.mass;
    this->time -= other.time;
    this->electric_current -= other.electric_current;
    this->thermodynamic_temperature -= other.thermodynamic_temperature;
    this->amount_of_substance -= other.amount_of_substance;
    this->luminous_intensity -= other.luminous_intensity;

    return *this;
}

bool
Dimensions
::operator!=(Dimensions const & other) const
{
    return !this->operator==(other);
}

Dimensions operator*(Dimensions l, Dimensions const & r)
{
    l *= r;
    return l;
}

Dimensions operator/(Dimensions l, Dimensions const & r)
{
    l /= r;
    return l;
}

std::ostream & operator<<(std::ostream & stream, Dimensions const & d)
{
    auto const print_dimension = [&](double dimension, std::string const & name) {
        if(dimension != 0)
        {
            stream << name;
            if(dimension != 1)
            {
                stream << "^" << dimension;
            }
            stream << " ";
        }
    };
    stream << "[ ";
    print_dimension(d.length, "L");
    print_dimension(d.mass, "M");
    print_dimension(d.time, "T");
    print_dimension(d.electric_current, "I");
    print_dimension(d.thermodynamic_temperature, "Θ");
    print_dimension(d.amount_of_substance, "N");
    print_dimension(d.luminous_intensity, "J");
    stream << "]";
    return stream;
}

Dimensions const Length{1, 0, 0, 0, 0, 0, 0};
Dimensions const Mass{0, 1, 0, 0, 0, 0, 0};
Dimensions const Time{0, 0, 1, 0, 0, 0, 0};
Dimensions const ElectricCurrent{0, 0, 0, 1, 0, 0, 0};
Dimensions const ThermodynamicTemperature{0, 0, 0, 0, 1, 0, 0};
Dimensions const AmountOfSubstance{0, 0, 0, 0, 0, 1, 0};
Dimensions const LuminousIntensity{0, 0, 0, 0, 0, 0, 1};

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

}

namespace std
{

sycomore::Dimensions pow(sycomore::Dimensions const & d, double s)
{
    return {
        d.length*s, d.mass*s, d.time*s, d.electric_current*s,
        d.thermodynamic_temperature*s, d.amount_of_substance*s,
        d.luminous_intensity*s};
}

}
