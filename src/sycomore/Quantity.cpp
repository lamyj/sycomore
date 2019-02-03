#include "Quantity.h"

#include <ostream>
#include <sstream>
#include <stdexcept>
#include "sycomore/Dimensions.h"

namespace sycomore
{

Quantity
::Quantity(double magnitude, Dimensions const & dimensions)
: magnitude(magnitude), dimensions(dimensions)
{
    // Nothing else.
}

bool
Quantity
::operator==(Quantity const & other) const
{
    return (
        this->magnitude == other.magnitude
        && this->dimensions == other.dimensions);
}

bool
Quantity
::operator!=(Quantity const & other) const
{
    return !this->operator==(other);
}

Quantity &
Quantity
::operator+=(Quantity const & other)
{
    if(this->dimensions != other.dimensions)
    {
        std::ostringstream message;
        message
            << "Unequal dimensions: "
            << this->dimensions << " != " << other.dimensions;
        throw std::runtime_error(message.str());
    }
    this->magnitude += other.magnitude;
    return *this;
}

Quantity &
Quantity
::operator-=(Quantity const & other)
{
    if(this->dimensions != other.dimensions)
    {
        std::ostringstream message;
        message
            << "Unequal dimensions: "
            << this->dimensions << " != " << other.dimensions;
        throw std::runtime_error(message.str());
    }
    this->magnitude -= other.magnitude;
    return *this;
}

Quantity &
Quantity
::operator*=(Quantity const & other)
{
    this->dimensions *= other.dimensions;
    this->magnitude *= other.magnitude;
    return *this;
}

Quantity &
Quantity
::operator*=(double scalar)
{
    this->magnitude *= scalar;
    return *this;
}

Quantity &
Quantity
::operator/=(Quantity const & other)
{
    this->dimensions /= other.dimensions;
    this->magnitude /= other.magnitude;
    return *this;
}

Quantity &
Quantity
::operator/=(double scalar)
{
    this->magnitude /= scalar;
    return *this;
}

double
Quantity
::convert_to(Quantity const & destination) const
{
    if(this->dimensions != destination.dimensions)
    {
        std::ostringstream message;
        message
            << "Unequal dimensions: "
            << this->dimensions << " != " << destination.dimensions;
        throw std::runtime_error(message.str());
    }
    return this->magnitude / destination.magnitude;
}

Quantity operator+(Quantity q) { return q; }
Quantity operator-(Quantity q) { q *= -1; return q; }

Quantity operator+(Quantity l, Quantity const & r) { l += r; return l; }
Quantity operator-(Quantity l, Quantity const & r) { l -= r; return l; }
Quantity operator*(Quantity l, Quantity const & r) { l *= r; return l; }
Quantity operator*(Quantity q, double s) { q *= s; return q; }
Quantity operator*(double s, Quantity q) { q *= s; return q; }
Quantity operator/(Quantity l, Quantity const & r) { l /= r; return l; }
Quantity operator/(Quantity q, double s) { q /= s; return q; }
Quantity operator/(double s, Quantity const & q)
{
    return Quantity{
        s/q.magnitude, {
            -q.dimensions.length,
            -q.dimensions.mass,
            -q.dimensions.time,
            -q.dimensions.electric_current,
            -q.dimensions.thermodynamic_temperature,
            -q.dimensions.amount_of_substance,
            -q.dimensions.luminous_intensity}};
}

std::ostream & operator<<(std::ostream & stream, Quantity const & q)
{
    stream << q.magnitude << " " << q.dimensions;
    return stream;
}

}
