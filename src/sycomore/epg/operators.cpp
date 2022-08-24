#include "operators.h"

#include <array>
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
pulse_magnetization_transfer(Real angle_a, Real phase_a, Real saturation)
{
    using std::cos; using std::exp; using std::pow;
    constexpr Complex const i{0,1};

    return {
        PULSE_MATRIX(angle_a, phase_a),
        exp(-saturation)
    };
}

std::pair<Real, Real> 
relaxation_single_pool(Real R1, Real R2, Real duration)
{
    auto const E_1 = std::exp(-duration*R1);
    auto const E_2 = std::exp(-duration*R2);
    return std::make_pair(E_1, E_2);
}

std::tuple<std::array<Complex, 8>, std::array<Real, 4>, std::array<Real, 2>>
relaxation_exchange(
    Real R1_a, Real R2_a, Real R1_b, Real R2_b,
    Real k_a, Real k_b, Real delta_b,
    Real M0_a, Real M0_b,
    Real duration)
{
    auto const exp_Lambda_T = relaxation_and_exchange_T(
        R2_a, R2_b, k_a, k_b, delta_b, duration);
    auto const exp_Lambda_L = relaxation_and_exchange_L(R1_a, R1_b, k_a, k_b, duration);
    auto const recovery = relaxation_and_exchange_recovery(
        R1_a, R1_b, k_a, k_b, M0_a, M0_b, exp_Lambda_L);
    return {exp_Lambda_T, exp_Lambda_L, recovery};
}

std::tuple<Real, std::array<Real, 4>, std::array<Real, 2>>
relaxation_magnetization_transfer(
    Real R1_a, Real R2_a, Real R1_b,
    Real k_a, Real k_b,
    Real M0_a, Real M0_b,
    Real duration)
{
    auto const exp_Lambda_T = std::exp(-R2_a*duration);
    auto const exp_Lambda_L = relaxation_and_exchange_L(
        R1_a, R1_b, k_a, k_b, duration);
    auto const recovery = relaxation_and_exchange_recovery(
        R1_a, R1_b, k_a, k_b, M0_a, M0_b, exp_Lambda_L);
    return {exp_Lambda_T, exp_Lambda_L, recovery};
}

std::array<Complex, 8>
relaxation_and_exchange_T(
    Real R2_a, Real R2_b, Real k_a, Real k_b, Real delta_b, Real duration)
{
    constexpr Complex const i{0,1};
    
    Real const a = duration*(-R2_a - k_a);
    Real const b = duration*(k_b);
    Real const c = duration*(k_a);
    Complex const d = duration*(-R2_b - k_b - 2*M_PI*i*delta_b);
    Complex const e = duration*(-R2_b - k_b + 2*M_PI*i*delta_b);
    
    if(b == 0. && c == 0.)
    {
        return {std::exp(a), 0, std::exp(a), 0, 0, std::exp(d), 0, std::exp(e)};
    }
    else if(b != 0. && c == 0.)
    {
        if(a == d && a == e)
        {
            auto const e_a = std::exp(a);
            return {e_a, b*e_a, e_a, b*e_a, 0, e_a, 0, e_a};
        }
        else if(a == d)
        {
            auto const e_a = std::exp(a);
            auto const e_e = std::exp(e);
            return {e_a, b*e_a, e_a, b*(e_a-e_e)/(a-e), 0, e_a, 0, e_e};
        }
        else if(a == e)
        {
            auto const e_a = std::exp(a);
            auto const e_d = std::exp(d);
            return {e_a, b*(e_a-e_d)/(a-d), e_a, b*e_a, 0, e_d, 0, e_a};
        }
        else
        {
            auto const e_a = std::exp(a);
            auto const e_d = std::exp(d);
            auto const e_e = std::exp(e);
            
            return {
                e_a, b*(e_a-e_d)/(a-d), e_a, b*(e_a-e_e)/(a-e), 0, e_d, 0, e_e};
        }
    }
    else
    {
        auto const Delta_1 = std::pow(a-d, 2)+4*b*c;
        auto const Delta_2 = std::pow(a-e, 2)+4*b*c;
        
        auto const l1 = (-std::sqrt(Delta_1) + a+d)/2.;
        auto const l2 = (+std::sqrt(Delta_1) + a+d)/2.;
        auto const l3 = (-std::sqrt(Delta_2) + a+e)/2.;
        auto const l4 = (+std::sqrt(Delta_2) + a+e)/2.;
        
        auto const e1 = std::exp(l1);
        auto const e2 = std::exp(l2);
        auto const e3 = std::exp(l3);
        auto const e4 = std::exp(l4);
        
        auto const dl_12 = (l1-l2);
        auto const dl_34 = (l3-l4);
        
        return {
            ((-d+l1)*e1+(d-l2)*e2)/(l1-l2), ((d-l1)*(d-l2)*(-e1+e2))/(c*(l1-l2)),
            ((-e+l3)*e3+(e-l4)*e4)/(l3-l4), ((e-l3)*(e-l4)*(-e3+e4))/(c*(l3-l4)),
            c*(e1-e2)/(l1-l2), ((-d+l1)*e2+(d-l2)*e1)/(l1-l2),
            c*(e3-e4)/(l3-l4), ((-e+l3)*e4+(e-l4)*e3)/(l3-l4)};
    }
}

std::array<Real, 4>
relaxation_and_exchange_L(
    Real R1_a, Real R1_b, Real k_a, Real k_b, Real duration)
{
    std::array<Real, 4> const A{
        duration*(-R1_a-k_a),       duration*(k_b),
              duration*(k_a), duration*(-R1_b-k_b)};
    return expm(A);
}

std::array<Real, 2>
relaxation_and_exchange_recovery(
    Real R1_a, Real R1_b, Real k_a, Real k_b, Real M0_a, Real M0_b,
    std::array<Real, 4> const & Xi_L)
{
    std::array<Real, 4> const Xi_L_minus_I{
        Xi_L[0]-1, Xi_L[1],
        Xi_L[2]  , Xi_L[3]-1};
        
    auto const det_Lambda_L = R1_a*R1_b + R1_a*k_b + R1_b*k_a;
    std::array<Real, 4> const inv_Lambda_L{
        (-R1_b-k_b)/det_Lambda_L,        -k_b/det_Lambda_L,
               -k_a/det_Lambda_L, (-R1_a-k_a)/det_Lambda_L};
    
    std::array<Real, 4> const product{
        Xi_L_minus_I[0]*inv_Lambda_L[0]+Xi_L_minus_I[1]*inv_Lambda_L[2],
        Xi_L_minus_I[0]*inv_Lambda_L[1]+Xi_L_minus_I[1]*inv_Lambda_L[3],
        Xi_L_minus_I[2]*inv_Lambda_L[0]+Xi_L_minus_I[3]*inv_Lambda_L[2],
        Xi_L_minus_I[2]*inv_Lambda_L[1]+Xi_L_minus_I[3]*inv_Lambda_L[3]};
    
    std::array<Real, 2> const C{R1_a*M0_a, R1_b*M0_b};
    
    return {
        product[0]*C[0]+product[1]*C[1],
        product[2]*C[0]+product[3]*C[1]};
}

std::array<Real, 4> expm(std::array<Real, 4> const & A)
{
    auto const Delta = std::pow(A[0]-A[3], 2) + 4*A[1]*A[2];
    
    static std::array<Real,4> const I{
        1, 0,
        0, 1};
    
    auto const lambda = (A[0]+A[3])/2.;
    
    std::array<Real, 4> result;
    
    if(Delta == 0.)
    {
        std::transform(
            A.begin(), A.end(), I.begin(), result.begin(),
            [&](Real a, Real i) {
                return std::exp(lambda) * (a + (1-lambda)*i); });
    }
    else if(Delta < 0)
    {
        auto const mu = std::sqrt(-Delta)/2.;
        std::transform(
            A.begin(), A.end(), I.begin(), result.begin(),
            [&](Real a, Real i) {
                return (
                    std::exp(lambda)
                    * (
                        (std::cos(mu)-lambda*std::sin(mu)/mu) * i
                        + std::sin(mu)/mu*a)); });
    }
    else
    {
        auto const nu = std::sqrt(Delta)/2.;
        
        if(std::fabs(lambda) < 1e-6)
        {
            std::transform(
                A.begin(), A.end(), I.begin(), result.begin(),
                [&](Real a, Real i) {
                    return std::cosh(nu)*i + std::sinh(nu)/nu*a; });
        }
        else if(lambda == nu || lambda == -nu)
        {
            std::transform(
                A.begin(), A.end(), I.begin(), result.begin(),
                [&](Real a, Real i) {
                    return i + (std::exp(2*lambda)-1.)/(2*lambda)*a; });
        }
        else
        {
            auto const xi1 = lambda + nu, xi2 = lambda - nu;
            std::transform(
                A.begin(), A.end(), I.begin(), result.begin(),
                [&](Real a, Real i) {
                    return (
                        (
                            (xi2*std::exp(xi1) - xi1*std::exp(xi2)) * i
                            + (std::exp(xi2)-std::exp(xi1)) * a
                        )
                        / (xi2-xi1)
                    ); });
        }
    }
    
    return result;
}

std::pair<Complex, Complex> phase_accumulation(Real angle)
{
    constexpr Complex const i{0,1};
    return {std::exp(i*angle), std::exp(i*-angle)};
}

}
    
}

}
