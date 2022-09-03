#define BOOST_TEST_MODULE isochromat_Model
#include <boost/test/unit_test.hpp>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#include "sycomore/isochromat/Model.h"

#include <xtensor/xio.hpp>

struct Fixture
{
    sycomore::isochromat::Positions positions;
    
    Fixture()
    : positions{xt::zeros<sycomore::isochromat::Position::value_type>({5, 3})}
    {
        for(std::size_t i=0; i<this->positions.shape()[0]; ++i)
        {
            auto const x = 2. * i / (this->positions.shape()[0]-1) - 1.;
            this->positions(i, 0) = x;
        }
    }
};

BOOST_AUTO_TEST_CASE(PulseUniformScratch)
{
    sycomore::isochromat::Positions positions{{0,0,0}};
    sycomore::isochromat::Model model(1, 0.1, {0,0,1}, positions);
    
    auto const op = model.build_pulse(M_PI/2, M_PI/3);
    sycomore::isochromat::Operator pulse{
        {{           0.25, std::sqrt(3)/4,  std::sqrt(3)/2, 0},
         { std::sqrt(3)/4,           0.75,            -0.5, 0},
         {-std::sqrt(3)/2,            0.5,               0, 0},
         {              0,              0,               0, 1}}
    };
    BOOST_TEST(xt::allclose(op, pulse));
}

BOOST_AUTO_TEST_CASE(PulseVariableScratch)
{
    sycomore::isochromat::Positions positions{{-1,0,0}, {1,0,0}};
    sycomore::isochromat::Model model(1, 0.1, {0,0,1}, positions);
    auto const op = model.build_pulse({M_PI/2., M_PI/3}, {M_PI/3., M_PI/2.});
    sycomore::isochromat::Operator pulse{
        {{           0.25, std::sqrt(3)/4, std::sqrt(3)/2, 0},
         { std::sqrt(3)/4,           0.75,           -0.5, 0},
         {-std::sqrt(3)/2,            0.5,              0, 0},
         {              0,              0,              0, 1}},
         
         {{           0.5, 0, std::sqrt(3)/2, 0},
         {              0, 1,              0, 0},
         {-std::sqrt(3)/2, 0,            0.5, 0},
         {              0, 0,              0, 1}}
    };
    BOOST_TEST(xt::allclose(op, pulse));
}

BOOST_AUTO_TEST_CASE(PulseUniformCombine)
{
    sycomore::isochromat::Positions positions{{0,0,0}};
    sycomore::isochromat::Model model(1, 0.1, {0,0,1}, positions);
    
    auto op = model.build_pulse(M_PI/2, M_PI/3);
    model.build_pulse(M_PI/3, M_PI/2, op);
    
    sycomore::isochromat::Operator pulse{
        {{           -0.625, 3*std::sqrt(3)/8, std::sqrt(3)/4, 0},
         {   std::sqrt(3)/4,             0.75,           -0.5, 0},
         {-3*std::sqrt(3)/8,           -0.125,          -0.75, 0},
         {                0,                0,              0, 1}}
    };
    BOOST_TEST(xt::allclose(op, pulse));
}

BOOST_AUTO_TEST_CASE(Relaxation)
{
    sycomore::isochromat::Positions positions{{{0,0,0}, {1,0,0}}};
    sycomore::isochromat::Model model(
        {1., 2.}, {0.1, 0.2}, {{0,0,2}, {0,0,1}}, positions);
    
    auto op = model.build_relaxation(1e-3);
    
    std::vector<sycomore::Real> E1{std::exp(-1e-3/1.), std::exp(-1e-3/2.)};
    std::vector<sycomore::Real> E2{std::exp(-1e-3/0.1), std::exp(-1e-3/0.2)};
    sycomore::isochromat::Operator relaxation{
        {{ E2[0],     0,     0,           0},
         {     0, E2[0],     0,           0},
         {     0,     0, E1[0], 2*(1-E1[0])},
         {     0,     0,     0,           1}},
        {{ E2[1],     0,     0,           0},
         {     0, E2[1],     0,           0},
         {     0,     0, E1[1], 1*(1-E1[1])},
         {     0,     0,     0,           1}},
    };
    BOOST_TEST(xt::allclose(op, relaxation));
}

BOOST_AUTO_TEST_CASE(PhaseAccumulation)
{
    sycomore::isochromat::Positions positions{{{0,0,0}, {1,0,0}}};
    sycomore::isochromat::Model model(1., 0.1, {0,0,1}, positions);
    
    auto op = model.build_phase_accumulation({M_PI/6., M_PI/3.});
    
    sycomore::isochromat::Operator phase_accumulation{
        {{ std::sqrt(3.)/2.,             -0.5, 0, 0},
         {              0.5, std::sqrt(3.)/2., 0, 0},
         {                0,                0, 1, 0},
         {                0,                0, 0, 1}},
        {{              0.5, -std::sqrt(3.)/2., 0, 0},
         { std::sqrt(3.)/2.,               0.5, 0, 0},
         {                0,                 0, 1, 0},
         {                0,                 0, 0, 1}},
    };
    BOOST_TEST(xt::allclose(op, phase_accumulation));
}
