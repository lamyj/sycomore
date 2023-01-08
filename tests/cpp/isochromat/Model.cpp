#define BOOST_TEST_MODULE isochromat_Model
#include <boost/test/unit_test.hpp>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#include "sycomore/isochromat/Model.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(PulseUniform)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model model(1*s, 0.1*s, {0,0,1}, {{0*m,0*m,0*m}});
    
    auto const op = model.build_pulse(M_PI/2*rad, M_PI/3*rad);
    sycomore::isochromat::Operator::Array pulse{
        {{           0.25, std::sqrt(3)/4,  std::sqrt(3)/2, 0},
         { std::sqrt(3)/4,           0.75,            -0.5, 0},
         {-std::sqrt(3)/2,            0.5,               0, 0},
         {              0,              0,               0, 1}}
    };
    BOOST_TEST(xt::allclose(op.array(), pulse));
}

BOOST_AUTO_TEST_CASE(PulseVariable)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model model(
        1*s, 0.1*s, {0,0,1}, {{-1*m,0*m,0*m}, {1*m,0*m,0*m}});
    auto const op = model.build_pulse(
        {M_PI/2.*rad, M_PI/3*rad}, {M_PI/3.*rad, M_PI/2.*rad});
    sycomore::isochromat::Operator::Array pulse{
        {{           0.25, std::sqrt(3)/4, std::sqrt(3)/2, 0},
         { std::sqrt(3)/4,           0.75,           -0.5, 0},
         {-std::sqrt(3)/2,            0.5,              0, 0},
         {              0,              0,              0, 1}},
         
         {{           0.5, 0, std::sqrt(3)/2, 0},
         {              0, 1,              0, 0},
         {-std::sqrt(3)/2, 0,            0.5, 0},
         {              0, 0,              0, 1}}
    };
    BOOST_TEST(xt::allclose(op.array(), pulse));
}

BOOST_AUTO_TEST_CASE(Relaxation)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model model(
        {1*s, 2*s}, {0.1*s, 0.2*s}, {{0,0,2}, {0,0,1}}, 
        {{0*m,0*m,0*m}, {1*m,0*m,0*m}});
    
    auto op = model.build_relaxation(1*ms);
    
    std::vector<sycomore::Real> E1{std::exp(-1e-3/1.), std::exp(-1e-3/2.)};
    std::vector<sycomore::Real> E2{std::exp(-1e-3/0.1), std::exp(-1e-3/0.2)};
    sycomore::isochromat::Operator::Array relaxation{
        {{ E2[0],     0,     0,           0},
         {     0, E2[0],     0,           0},
         {     0,     0, E1[0], 2*(1-E1[0])},
         {     0,     0,     0,           1}},
        {{ E2[1],     0,     0,           0},
         {     0, E2[1],     0,           0},
         {     0,     0, E1[1], 1*(1-E1[1])},
         {     0,     0,     0,           1}},
    };
    BOOST_TEST(xt::allclose(op.array(), relaxation));
}

BOOST_AUTO_TEST_CASE(PhaseAccumulation)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model model(
        1*s, 0.1*s, {0,0,1}, {{0*m,0*m,0*m}, {1*m,0*m,0*m}});
    
    auto op = model.build_phase_accumulation({M_PI/6.*rad, M_PI/3.*rad});
    
    sycomore::isochromat::Operator::Array phase_accumulation{
        {{ std::sqrt(3.)/2.,             -0.5, 0, 0},
         {              0.5, std::sqrt(3.)/2., 0, 0},
         {                0,                0, 1, 0},
         {                0,                0, 0, 1}},
        {{              0.5, -std::sqrt(3.)/2., 0, 0},
         { std::sqrt(3.)/2.,               0.5, 0, 0},
         {                0,                 0, 1, 0},
         {                0,                 0, 0, 1}},
    };
    BOOST_TEST(xt::allclose(op.array(), phase_accumulation));
}

BOOST_AUTO_TEST_CASE(TimeIntervalUniform)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model model(
        1*s, 0.1*s, {0, 0, 1}, {{0*m, 0*m, 0*m}, {1*mm, 2*mm, 3*mm}});
    
    auto op = model.build_time_interval(
        10_ms, 400*Hz, {20*mT/m, 0*mT/m, 10*mT/m});
    
    auto combined = model.build_phase_accumulation(2*M_PI*400*10e-3);
    combined.preMultiply(model.build_relaxation(10*ms));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 0), combined.array()));
    
    combined = model.build_phase_accumulation(
        2*M_PI*400*10e-3 + sycomore::gamma.magnitude*50e-6*10e-3);
    combined.preMultiply(model.build_relaxation(10*ms));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 1), combined.array()));
}

BOOST_AUTO_TEST_CASE(TimeIntervalVariable)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model model(
        1*s, 0.1*s, {0, 0, 1}, {{-1*mm, 2*mm, 3*mm}, {1*mm, 0*mm, 3*mm}});
    
    auto op = model.build_time_interval(
        10*ms, {400*Hz, 600*Hz},
        {{20*mT/m, 0*mT/m, 10*mT/m}, {15*mT/m, 17e-3*mT/m, 0*mT/m}});
    
    auto combined = model.build_phase_accumulation(
        2*M_PI*400*10e-3 + sycomore::gamma.magnitude*10e-6*10e-3);
    combined.preMultiply(model.build_relaxation(10*ms));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 0), combined.array()));
    
    combined = model.build_phase_accumulation(
        2*M_PI*600*10e-3 + sycomore::gamma.magnitude*15e-6*10e-3);
    combined.preMultiply(model.build_relaxation(10*ms));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 1), combined.array()));
}

BOOST_AUTO_TEST_CASE(T1)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model const model_1(
        1*ms, 2*ms, {3., 4., 5.}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((model_1.T1() == sycomore::TensorQ<1>{{1*ms, 1*ms}}));
    
    sycomore::isochromat::Model const model_2(
        {1*ms, 2*ms}, {3*ms, 4*ms},
        {{5., 6., 7.}, {8., 9., 10.}}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((model_2.T1() == sycomore::TensorQ<1>{{1*ms, 2*ms}}));
}

BOOST_AUTO_TEST_CASE(T2)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model const model_1(
        1*ms, 2*ms, {3., 4., 5.}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((model_1.T2() == sycomore::TensorQ<1>{{2*ms, 2*ms}}));
    
    sycomore::isochromat::Model const model_2(
        {1*ms, 2*ms}, {3*ms, 4*ms},
        {{5., 6., 7.}, {8., 9., 10.}}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((model_2.T2() == sycomore::TensorQ<1>{{3*ms, 4*ms}}));
}

BOOST_AUTO_TEST_CASE(M0)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model const model_1(
        1*ms, 2*ms, {3., 4., 5.}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((model_1.M0() == sycomore::TensorR<1>{{5., 5.}}));
    
    sycomore::isochromat::Model const model_2(
        {1*ms, 2*ms}, {3*ms, 4*ms},
        {{5., 6., 7.}, {8., 9., 10.}}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((model_2.M0() == sycomore::TensorR<1>{{7., 10.}}));
}

BOOST_AUTO_TEST_CASE(Magnetization)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model const model_1(
        1*s, 2*s, {3., 4., 5.}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((
        model_1.magnetization()
        == sycomore::TensorR<2>{{3., 4., 5.}, {3., 4., 5.}}));
    
    sycomore::isochromat::Model const model_2(
        {1*s, 2*s}, {3*s, 4*s},
        {{5., 6., 7.}, {8., 9., 10.}}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((
        model_2.magnetization()
        == sycomore::TensorR<2>{{5., 6., 7.}, {8., 9., 10.}}));
}

BOOST_AUTO_TEST_CASE(Positions)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Model const model(
        1*ms, 2*ms, {3., 4., 5.}, {{0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}});
    BOOST_TEST((
        model.positions() == sycomore::TensorQ<2>{
            {0*m, 0*m, 11*m}, {0*m, 0*m, 12*m}}));
}

BOOST_AUTO_TEST_CASE(Apply)
{
    using namespace sycomore::units;
    
    sycomore::isochromat::Operator operator_({
        {
            {1, 2, 3, 10},
            {4, 5, 6, 20},
            {7, 8, 9, 30},
            {0, 0, 0, 1},
        },
        {
            {10, 11, 12, 40},
            {13, 14, 15, 50},
            {16, 17, 18, 60},
            {0, 0, 0, 1},
        }
    });
    
    sycomore::isochromat::Model model(
        {1*s, 1*s}, {1*s, 1*s},
        {{19., 20., 21.}, {22., 23., 24.}}, {{0*m, 0*m, 0*m}, {0*m, 0*m, 1*m}});
    model.apply(operator_);
    
    decltype(model.magnetization()) magnetization{
        {132, 322, 512}, {801, 1018, 1235}};
    
    BOOST_TEST(xt::allclose(model.magnetization(), magnetization));
}
