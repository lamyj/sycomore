#include "sycomore.h"

namespace sycomore
{

Real rad2deg(Real const value)
{
    return value/M_PI*180;
}

Real deg2rad(Real const value)
{
    return value/180*M_PI;
}

}
