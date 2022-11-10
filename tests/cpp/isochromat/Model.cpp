#define BOOST_TEST_MODULE isochromat_Model
#include <boost/test/unit_test.hpp>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#include "sycomore/isochromat/Model.h"

BOOST_AUTO_TEST_CASE(PulseUniform)
{
    xt::xtensor<sycomore::Real, 2> positions{{0,0,0}};
    sycomore::isochromat::Model model(1, 0.1, {0,0,1}, positions);
    
    auto const op = model.build_pulse(M_PI/2, M_PI/3);
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
    xt::xtensor<sycomore::Real, 2> positions{{-1,0,0}, {1,0,0}};
    sycomore::isochromat::Model model(1, 0.1, {0,0,1}, positions);
    auto const op = model.build_pulse({M_PI/2., M_PI/3}, {M_PI/3., M_PI/2.});
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
    xt::xtensor<sycomore::Real, 2> positions{{{0,0,0}, {1,0,0}}};
    sycomore::isochromat::Model model(
        {1., 2.}, {0.1, 0.2}, {{0,0,2}, {0,0,1}}, positions);
    
    auto op = model.build_relaxation(1e-3);
    
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
    xt::xtensor<sycomore::Real, 2> positions{{{0,0,0}, {1,0,0}}};
    sycomore::isochromat::Model model(1., 0.1, {0,0,1}, positions);
    
    auto op = model.build_phase_accumulation({M_PI/6., M_PI/3.});
    
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
    xt::xtensor<sycomore::Real, 2> positions{
        {0, 0, 0}, {1e-3, 2e-3, 3e-3}};
    sycomore::isochromat::Model model(1., 0.1, {0, 0, 1}, positions);
    
    auto op = model.build_time_interval(10e-3, 400, {20e-3, 0, 10e-3});
    
    auto combined = model.build_phase_accumulation(2*M_PI*400*10e-3);
    combined.preMultiply(model.build_relaxation(10e-3));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 0), combined.array()));
    
    combined = model.build_phase_accumulation(
        2*M_PI*400*10e-3 + sycomore::gamma.magnitude*50e-6*10e-3);
    combined.preMultiply(model.build_relaxation(10e-3));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 1), combined.array()));
}

BOOST_AUTO_TEST_CASE(TimeIntervalVariable)
{
    xt::xtensor<sycomore::Real, 2> positions{
        {-1e-3, 2e-3, 3e-3}, {1e-3, 0, 3e-3}};
    sycomore::isochromat::Model model(1., 0.1, {0, 0, 1}, positions);
    
    auto op = model.build_time_interval(
        10e-3, {400, 600}, {{20e-3, 0, 10e-3}, {15e-3, 17e-3, 0}});
    
    auto combined = model.build_phase_accumulation(
        2*M_PI*400*10e-3 + sycomore::gamma.magnitude*10e-6*10e-3);
    combined.preMultiply(model.build_relaxation(10e-3));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 0), combined.array()));
    
    combined = model.build_phase_accumulation(
        2*M_PI*600*10e-3 + sycomore::gamma.magnitude*15e-6*10e-3);
    combined.preMultiply(model.build_relaxation(10e-3));
    BOOST_TEST(xt::allclose(xt::view(op.array(), 1), combined.array()));
}

BOOST_AUTO_TEST_CASE(T1)
{
    sycomore::isochromat::Model const model_1(
        1., 2., {3., 4., 5.}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((model_1.T1() == xt::xtensor<double, 1>{{1.}}));
    
    sycomore::isochromat::Model const model_2(
        {1., 2.}, {3., 4.},
        {{5., 6., 7., 1.}, {8., 9., 10., 1.}}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((model_2.T1() == xt::xtensor<double, 1>{{1., 2.}}));
}

BOOST_AUTO_TEST_CASE(T2)
{
    sycomore::isochromat::Model const model_1(
        1., 2., {3., 4., 5.}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((model_1.T2() == xt::xtensor<double, 1>{{2.}}));
    
    sycomore::isochromat::Model const model_2(
        {1., 2.}, {3., 4.},
        {{5., 6., 7., 1.}, {8., 9., 10., 1.}}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((model_2.T2() == xt::xtensor<double, 1>{{3., 4.}}));
}

BOOST_AUTO_TEST_CASE(M0)
{
    sycomore::isochromat::Model const model_1(
        1., 2., {3., 4., 5.}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((model_1.M0() == xt::xtensor<double, 1>{{5.}}));
    
    sycomore::isochromat::Model const model_2(
        {1., 2.}, {3., 4.},
        {{5., 6., 7., 1.}, {8., 9., 10., 1.}}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((model_2.M0() == xt::xtensor<double, 1>{{7., 10.}}));
}

BOOST_AUTO_TEST_CASE(Magnetization)
{
    sycomore::isochromat::Model const model_1(
        1., 2., {3., 4., 5.}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((
        model_1.magnetization()
        == xt::xtensor<double, 2>{{3., 4., 5., 1.}, {3., 4., 5., 1.}}));
    
    sycomore::isochromat::Model const model_2(
        {1., 2.}, {3., 4.},
        {{5., 6., 7., 1.}, {8., 9., 10., 1.}}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((
        model_2.magnetization()
        == xt::xtensor<double, 2>{{5., 6., 7., 1.}, {8., 9., 10., 1.}}));
}

BOOST_AUTO_TEST_CASE(Positions)
{
    sycomore::isochromat::Model const model(
        1., 2., {3., 4., 5.}, {{0., 0., 11.}, {0., 0., 12.}});
    BOOST_TEST((
        model.positions()
        == xt::xtensor<double, 2>{{0., 0., 11.}, {0., 0., 12.}}));
}

BOOST_AUTO_TEST_CASE(Apply)
{
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
        {1., 1.}, {1., 1.},
        {{19., 20., 21., 1.}, {22., 23., 24., 1.}}, {{0., 0., 0.}, {0., 0., 1.}});
    model.apply(operator_);
    
    decltype(model.magnetization()) magnetization{
        {132, 322, 512, 1},
        {801, 1018, 1235, 1}};
    
    BOOST_TEST(xt::allclose(model.magnetization(), magnetization));
}
