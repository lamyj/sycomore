#include "units.h"

namespace sycomore
{

namespace units
{

template<int L, int M, int T, int I, int Theta, int N, int J>
Unit<L,M,T,I,Theta,N,J>
::Unit()
: Unit(0.0)
{
    // Nothing else.
}

template<int L, int M, int T, int I, int Theta, int N, int J>
Unit<L,M,T,I,Theta,N,J>
::Unit(double value)
: value(value)
{
    // Nothing else.
}

template<int L, int M, int T, int I, int Theta, int N, int J>
double
Unit<L,M,T,I,Theta,N,J>
::convert_to(Unit<L,M,T,I,Theta,N,J> dest) const
{
    return this->value / dest.value;
}

template<int ... Args>
Unit<Args ...>
operator+(Unit<Args ...> const & x, Unit<Args ...> const & y)
{
    return Unit<Args ...>(x.value+y.value);
}

template<int ... Args>
Unit<Args ...>
operator-(Unit<Args ...> const & x, Unit<Args ...> const & y)
{
    return Unit<Args ...>(x.value-y.value);
}

template<int ... Args>
Unit<Args ...>
operator*(double const & scalar, Unit<Args ...> const & unit)
{
    return Unit<Args ...>(scalar*unit.value);
}

template<int ... Args>
Unit<Args ...>
operator*(Unit<Args ...> const & unit, double const & scalar)
{
    return scalar*unit;
}

template<int ... Args>
pow<Unit<Args ...>, -1>
operator/(double const & scalar, Unit<Args ...> const & unit)
{
    return pow<Unit<Args ...>, -1>(scalar/unit.value);
}

template<int ... Args>
Unit<Args ...>
operator/(Unit<Args ...> const & unit, double const & scalar)
{
    return Unit<Args ...>(unit.value/scalar);
}

template<typename Unit1, typename Unit2>
mult<Unit1, Unit2>
operator*(Unit1 const & unit_1, Unit2 const & unit_2)
{
    return mult<Unit1, Unit2>(unit_1.value*unit_2.value);
}

template<typename Unit1, typename Unit2>
div<Unit1, Unit2>
operator/(Unit1 const & unit_1, Unit2 const & unit_2)
{
    return div<Unit1, Unit2>(unit_1.value/unit_2.value);
}

}

}
