#ifndef _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f
#define _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f

#include <ostream>
#include "sycomore/Dimensions.h"

namespace sycomore
{

/// @brief Quantity in the SI system.
class Quantity
{
public:
    /// @brief Magnitude of the quantity in base units.
    double magnitude;
    Dimensions dimensions;

    Quantity(double magnitude={}, Dimensions const & dimensions={});

    bool operator==(Quantity const & other) const;
    bool operator!=(Quantity const & other) const;

    Quantity & operator+=(Quantity const & other);
    Quantity & operator-=(Quantity const & other);
    Quantity & operator*=(Quantity const & other);
    Quantity & operator*=(double scalar);
    Quantity & operator/=(Quantity const & other);
    Quantity & operator/=(double scalar);
    Quantity & operator%=(Quantity const & other);
    Quantity & operator%=(double scalar);

    double convert_to(Quantity const & destination) const;
};

bool operator<(Quantity const & l, Quantity const & r);
bool operator<=(Quantity const & l, Quantity const & r);
bool operator>(Quantity const & l, Quantity const & r);
bool operator>=(Quantity const & l, Quantity const & r);

Quantity operator+(Quantity q);
Quantity operator-(Quantity q);
Quantity operator+(Quantity l, Quantity const & r);
Quantity operator-(Quantity l, Quantity const & r);
Quantity operator*(Quantity l, Quantity const & r);
Quantity operator*(Quantity q, double s);
Quantity operator*(double s, Quantity q);
Quantity operator/(Quantity l, Quantity const & r);
Quantity operator/(Quantity q, double s);
Quantity operator/(double s, Quantity const & q);
Quantity operator%(Quantity l, Quantity const & r);
Quantity operator%(Quantity q, double s);

std::ostream & operator<<(std::ostream & stream, Quantity const & q);

}

namespace std
{

sycomore::Quantity pow(sycomore::Quantity q, double e);

template<>
struct hash<sycomore::Quantity>
{
    size_t operator()(sycomore::Quantity const & q) const;
};

}

#endif // _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f
