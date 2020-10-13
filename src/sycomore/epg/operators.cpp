#include "operators.h"

#include <cmath>
#include <tuple>
#include <utility>
#include <vector>

#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

namespace epg
{

namespace operators
{

std::vector<Complex> pulse(Real angle, Real phase)
{
    using std::cos; using std::exp; using std::pow;
    
    auto const & a = angle;
    auto const & p = phase;
    
    constexpr Complex const i{0,1};

    return {
        pow(cos(a/2), 2),              exp(2.*i*p)*pow(sin(a/2), 2), -i*exp( i*p)*sin(a),
        exp(-2.*i*p)*pow(sin(a/2), 2), pow(cos(a/2), 2),              i*exp(-i*p)*sin(a),
        -i/2.*exp(-i*p)*sin(a),        i/2.*exp(i*p)*sin(a),          cos(a)
    };
}

std::pair<Real, Real> 
relaxation(Real R1, Real R2, Real duration)
{
    auto const E_1 = std::exp(-duration*R1);
    auto const E_2 = std::exp(-duration*R2);
    return std::make_pair(E_1, E_2);
}

std::pair<Complex, Complex> phase_accumulation(Real angle)
{
    constexpr Complex const i{0,1};
    return {std::exp(i*angle), std::exp(i*-angle)};
}

}
    
}

}
