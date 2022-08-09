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

#define PULSE_MATRIX(a, p) \
    pow(cos(a/2), 2),              exp(2.*i*p)*pow(sin(a/2), 2), -i*exp( i*p)*sin(a), \
    exp(-2.*i*p)*pow(sin(a/2), 2), pow(cos(a/2), 2),              i*exp(-i*p)*sin(a), \
    -i/2.*exp(-i*p)*sin(a),        i/2.*exp(i*p)*sin(a),          cos(a)

std::vector<Complex> pulse_single_pool(Real angle, Real phase)
{
    using std::cos; using std::exp; using std::pow;
    constexpr Complex const i{0,1};

    return { PULSE_MATRIX(angle, phase) };
}

std::vector<Complex>
pulse_exchange(Real angle_a, Real phase_a, Real angle_b, Real phase_b)
{
    using std::cos; using std::exp; using std::pow;
    constexpr Complex const i{0,1};

    return {
        PULSE_MATRIX(angle_a, phase_a),
        PULSE_MATRIX(angle_b, phase_b)
    };
}

std::vector<Complex>
pulse_magnetization_transfer(
    Real angle_a, Real phase_a, Real saturation_rate, Real duration)
{
    using std::cos; using std::exp; using std::pow;
    constexpr Complex const i{0,1};

    return {
        PULSE_MATRIX(angle_a, phase_a),
        exp(-saturation_rate*duration)
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
