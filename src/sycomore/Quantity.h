#ifndef _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f
#define _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f

#include <ostream>
#include "sycomore/Dimensions.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

/// @brief Quantity in the SI system.
class SYCOMORE_API Quantity
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

SYCOMORE_API bool operator<(Quantity const & l, Quantity const & r);
SYCOMORE_API bool operator<=(Quantity const & l, Quantity const & r);
SYCOMORE_API bool operator>(Quantity const & l, Quantity const & r);
SYCOMORE_API bool operator>=(Quantity const & l, Quantity const & r);

SYCOMORE_API Quantity operator+(Quantity q);
SYCOMORE_API Quantity operator-(Quantity q);
SYCOMORE_API Quantity operator+(Quantity l, Quantity const & r);
SYCOMORE_API Quantity operator-(Quantity l, Quantity const & r);
SYCOMORE_API Quantity operator*(Quantity l, Quantity const & r);
SYCOMORE_API Quantity operator*(Quantity q, double s);
SYCOMORE_API Quantity operator*(double s, Quantity q);
SYCOMORE_API Quantity operator/(Quantity l, Quantity const & r);
SYCOMORE_API Quantity operator/(Quantity q, double s);
SYCOMORE_API Quantity operator/(double s, Quantity const & q);
SYCOMORE_API Quantity operator%(Quantity l, Quantity const & r);
SYCOMORE_API Quantity operator%(Quantity q, double s);

SYCOMORE_API std::ostream & operator<<(std::ostream & stream, Quantity const & q);

}

namespace std
{

SYCOMORE_API sycomore::Quantity pow(sycomore::Quantity q, double e);

template<>
struct SYCOMORE_API hash<sycomore::Quantity>
{
    std::size_t operator()(sycomore::Quantity const & q) const;
};

}

#endif // _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f
