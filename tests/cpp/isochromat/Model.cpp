#define BOOST_TEST_MODULE isochromat_Model
#include <boost/test/unit_test.hpp>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#include "sycomore/isochromat/Model.h"

#include <xtensor/xio.hpp>

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
