#include "Dimensions.h"

#include <ostream>

namespace sycomore
{

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
    auto const print_dimension = [&](int dimension, std::string const & name) {
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

}

namespace std
{

sycomore::Dimensions pow(sycomore::Dimensions const & d, int s)
{
    return {
        d.length*s, d.mass*s, d.time*s, d.electric_current*s,
        d.thermodynamic_temperature*s, d.amount_of_substance*s,
        d.luminous_intensity*s};
}

}