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
    Quantity & operator+=(double s);
    
    Quantity & operator-=(Quantity const & other);
    Quantity & operator-=(double s);
    
    Quantity & operator*=(Quantity const & other);
    Quantity & operator*=(double scalar);
    Quantity & operator/=(Quantity const & other);
    Quantity & operator/=(double scalar);
    Quantity & operator%=(Quantity const & other);
    Quantity & operator%=(double scalar);

    /**
     * @brief Return the scalar value of the quantity converted to the given 
     * unit.
    */
    double convert_to(Quantity const & destination) const;
    
    operator double() const;
};

template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator==(Quantity const & q, T s) { return q == Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator==(T s, Quantity const & q) { return q == Quantity(s); }

template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator!=(Quantity const & q, T s) { return q != Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator!=(T s, Quantity const & q) { return q != Quantity(s); }

SYCOMORE_API bool operator<(Quantity const & l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<(Quantity const & q, T s) { return q < Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<(T s, Quantity const & q) { return Quantity(s) < q; }

SYCOMORE_API bool operator<=(Quantity const & l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<=(Quantity const & q, T s) { return q <= Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator<=(T s, Quantity const & q) { return Quantity(s) <= q; }

SYCOMORE_API bool operator>(Quantity const & l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>(Quantity const & q, T s) { return q > Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>(T s, Quantity const & q) { return Quantity(s) > q; }

SYCOMORE_API bool operator>=(Quantity const & l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>=(Quantity const & q, T s) { return q >= Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
bool operator>=(T s, Quantity const & q) { return Quantity(s) >= q; }

SYCOMORE_API Quantity operator+(Quantity q);
SYCOMORE_API Quantity operator-(Quantity q);

SYCOMORE_API Quantity operator+(Quantity l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator+(Quantity const & q, T s) { return q+Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator+(T s, Quantity const & q) { return Quantity(s)+q; }

SYCOMORE_API Quantity operator-(Quantity l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator-(Quantity const & q, T s) { return q-Quantity(s); }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator-(T s, Quantity const & q) { return Quantity(s)-q; }

SYCOMORE_API Quantity operator*(Quantity l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator*(Quantity q, T s) { q *= s; return q; }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator*(T s, Quantity const & q) { return q*s; }

SYCOMORE_API Quantity operator/(Quantity l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator/(Quantity q, T s) { q /= s; return q; }
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator/(T s, Quantity const & q) { return Quantity(s)/q; }

SYCOMORE_API Quantity operator%(Quantity l, Quantity const & r);
template<typename T, typename std::enable_if<std::is_arithmetic<T>::value, int>::type=0>
Quantity operator%(Quantity q, T s) { q %= double(s); return q; }

SYCOMORE_API std::ostream & operator<<(std::ostream & stream, Quantity const & q);

}

namespace std
{

SYCOMORE_API sycomore::Quantity abs(sycomore::Quantity q);
SYCOMORE_API sycomore::Quantity pow(sycomore::Quantity q, double e);

SYCOMORE_API sycomore::Quantity round(sycomore::Quantity q);
SYCOMORE_API sycomore::Quantity trunc(sycomore::Quantity q);
SYCOMORE_API sycomore::Quantity floor(sycomore::Quantity q);
SYCOMORE_API sycomore::Quantity ceil(sycomore::Quantity q);

template<>
struct SYCOMORE_API hash<sycomore::Quantity>
{
    std::size_t operator()(sycomore::Quantity const & q) const;
};

}

#endif // _bd3de17b_e4fa_4e7f_8d72_8ac9df01606f
