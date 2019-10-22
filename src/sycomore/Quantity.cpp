#include "Quantity.h"

#include <cmath>
#include <ostream>
#include <sstream>
#include <stdexcept>

#include "sycomore/Dimensions.h"
#include "sycomore/hash.h"

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
            << "Addition requires equal dimensions: "
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
            << "Subtraction requires equal dimensions: "
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

Quantity &
Quantity
::operator%=(Quantity const & other)
{
    if(this->dimensions != other.dimensions)
    {
        std::ostringstream message;
        message
            << "Modulo requires equal dimensions: "
            << this->dimensions << " != " << other.dimensions;
        throw std::runtime_error(message.str());
    }
    this->magnitude = std::fmod(this->magnitude, other.magnitude);
    return *this;
}

Quantity &
Quantity
::operator%=(double scalar)
{
    this->magnitude = std::fmod(this->magnitude, scalar);
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
            << "Conversion requires equal dimensions: "
            << this->dimensions << " != " << destination.dimensions;
        throw std::runtime_error(message.str());
    }
    return this->magnitude / destination.magnitude;
}

bool operator<(Quantity const & l, Quantity const & r)
{
    if(l.dimensions != r.dimensions)
    {
        std::ostringstream message;
        message
            << "Comparison requires equal dimensions: "
            << l.dimensions << " != " << r.dimensions;
        throw std::runtime_error(message.str());
    }
    return l.magnitude < r.magnitude;
}

bool operator<=(Quantity const & l, Quantity const & r)
{
    return l<r || l==r;
}

bool operator>(Quantity const & l, Quantity const & r)
{
    return !(l<=r);
}

bool operator>=(Quantity const & l, Quantity const & r)
{
    return !(l<r);
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
Quantity operator%(Quantity q, double s) { q %= s; return q; }
Quantity operator%(Quantity l, Quantity const & r) { l %= r; return l; }

std::ostream & operator<<(std::ostream & stream, Quantity const & q)
{
    stream << q.magnitude << " " << q.dimensions;
    return stream;
}

}

namespace std
{

sycomore::Quantity pow(sycomore::Quantity q, double e)
{
    q.dimensions = std::pow(q.dimensions, e);
    q.magnitude = std::pow(q.magnitude, e);
    return q;
}

//template<>
std::size_t
hash<sycomore::Quantity>
::operator()(sycomore::Quantity const & q) const
{
    std::size_t seed=0;
    hash<double> hasher;
    sycomore::combine_hashes(seed, hasher(q.magnitude));
    sycomore::combine_hashes(seed, hasher(q.dimensions.length));
    sycomore::combine_hashes(seed, hasher(q.dimensions.mass));
    sycomore::combine_hashes(seed, hasher(q.dimensions.time));
    sycomore::combine_hashes(seed, hasher(q.dimensions.electric_current));
    sycomore::combine_hashes(seed, hasher(q.dimensions.thermodynamic_temperature));
    sycomore::combine_hashes(seed, hasher(q.dimensions.amount_of_substance));
    sycomore::combine_hashes(seed, hasher(q.dimensions.luminous_intensity));

    return seed;
}

}
