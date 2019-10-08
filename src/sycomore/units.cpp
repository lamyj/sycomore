#include "units.h"

#include <cmath>
#include <ostream>

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"

// WARNING 3.6.2-2: the order of initialization in different translation units 
// is undefined. Avoid using pre-defined dimensions.

namespace sycomore
{

namespace units
{

#define SYCOMORE_DEFINE_UNIT(dimensions, name, factor) \
    Quantity const name{factor, dimensions}; \
    Quantity operator "" _##name(unsigned long long v) { return double(v)*name; } \
    Quantity operator "" _##name(long double v) { return double(v)*name; }

#define SYCOMORE_DEFINE_UNITS(Type, name) \
    SYCOMORE_DEFINE_UNIT(Type, Y##name, 1e24) \
    SYCOMORE_DEFINE_UNIT(Type, Z##name, 1e21) \
    SYCOMORE_DEFINE_UNIT(Type, E##name, 1e18) \
    SYCOMORE_DEFINE_UNIT(Type, P##name, 1e15) \
    SYCOMORE_DEFINE_UNIT(Type, T##name, 1e12) \
    SYCOMORE_DEFINE_UNIT(Type, G##name, 1e9) \
    SYCOMORE_DEFINE_UNIT(Type, M##name, 1e6) \
    SYCOMORE_DEFINE_UNIT(Type, k##name, 1e3) \
    SYCOMORE_DEFINE_UNIT(Type, h##name, 1e2) \
    SYCOMORE_DEFINE_UNIT(Type, da##name, 1e1) \
    SYCOMORE_DEFINE_UNIT(Type, name, 1) \
    SYCOMORE_DEFINE_UNIT(Type, d##name, 1e-1) \
    SYCOMORE_DEFINE_UNIT(Type, c##name, 1e-2) \
    SYCOMORE_DEFINE_UNIT(Type, m##name, 1e-3) \
    SYCOMORE_DEFINE_UNIT(Type, u##name, 1e-6) \
    SYCOMORE_DEFINE_UNIT(Type, n##name, 1e-9) \
    SYCOMORE_DEFINE_UNIT(Type, p##name, 1e-12) \
    SYCOMORE_DEFINE_UNIT(Type, f##name, 1e-15) \
    SYCOMORE_DEFINE_UNIT(Type, a##name, 1e-18) \
    SYCOMORE_DEFINE_UNIT(Type, z##name, 1e-21) \
    SYCOMORE_DEFINE_UNIT(Type, y##name, 1e-24)

SYCOMORE_DEFINE_UNITS((Dimensions{1, 0, 0, 0, 0, 0, 0}), m)
// WARNING: the base unit is kg, but since we only define the cast operators,
// this simpler form works.
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), Yg, 1e21)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), Zg, 1e18)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), Eg, 1e15)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), Pg, 1e12)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), Tg, 1e9)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), Gg, 1e6)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), Mg, 1e3)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), kg, 1)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), hg, 1e-1)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), dag, 1e-2)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), g, 1e-3)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), dg, 1e-4)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), cg, 1e-5)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), mg, 1e-6)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), ug, 1e-9)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), ng, 1e-12)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), pg, 1e-15)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), fg, 1e-18)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), ag, 1e-21)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), zg, 1e-24)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, 0, 0, 0, 0, 0}), yg, 1e-27)

SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 1, 0, 0, 0, 0}), s)
//SYCOMORE_DEFINE_UNIT(Time, min)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 0, 1, 0, 0, 0, 0}), h, 3600)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 0, 1, 0, 0, 0}), A)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 0, 0, 1, 0, 0}), K)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 0, 0, 0, 1, 0}), mol)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 0, 0, 0, 0, 1}), cd)

SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 0, 0, 0, 0, 0}), rad)
SYCOMORE_DEFINE_UNIT((Dimensions{0, 0, 0, 0, 0, 0, 0}), deg, M_PI/180.)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 0, 0, 0, 0, 0}), sr)

SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, -1, 0, 0, 0, 0}), Hz)
SYCOMORE_DEFINE_UNITS((Dimensions{1, 1, -2, 0, 0, 0, 0}), N)
SYCOMORE_DEFINE_UNITS((Dimensions{-1, 1, -2, 0, 0, 0, 0}), Pa)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 1, -2, 0, 0, 0, 0}), J)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 1, -3, 0, 0, 0, 0}), W)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 1, 1, 0, 0, 0}), C)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 1, -3, -1, 0, 0, 0}), V)
SYCOMORE_DEFINE_UNITS((Dimensions{-2, -1, 4, 2, 0, 0, 0}), F)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 1, -3, -2, 0, 0, 0}), Ohm)
SYCOMORE_DEFINE_UNITS((Dimensions{-2, -1, 3, 2, 0, 0, 0}), S)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 1, -2, -1, 0, 0, 0}), Wb)

// WARNING: _MT is a macro on Windows
#ifdef _WIN32
#pragma push_macro("_MT")
#undef _MT
#endif

SYCOMORE_DEFINE_UNITS((Dimensions{0, 1, -2, -1, 0, 0, 0}), T);

#ifdef _WIN32
#pragma pop_macro("_MT")
#endif

SYCOMORE_DEFINE_UNIT((Dimensions{0, 1, -2, -1, 0, 0, 0}), G, 1e-4)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 1, -2, -2, 0, 0, 0}), H)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, 0, 0, 0, 0, 1}), lm)
SYCOMORE_DEFINE_UNITS((Dimensions{-2, 0, 0, 0, 0, 0, 1}), lx)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, -1, 0, 0, 0, 0}), Bq)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 0, -2, 0, 0, 0, 0}), Gy)
SYCOMORE_DEFINE_UNITS((Dimensions{2, 0, -2, 0, 0, 0, 0}), Sv)
SYCOMORE_DEFINE_UNITS((Dimensions{0, 0, -1, 0, 0, 1, 0}), kat)

}

}
