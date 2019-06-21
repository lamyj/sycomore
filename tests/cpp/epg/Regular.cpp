#define BOOST_TEST_MODULE epg_Regular
#include <boost/test/unit_test.hpp>

#include "sycomore/epg/Regular.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

#define TEST_COMPLEX_EQUAL(v1, v2) \
    { \
        sycomore::Complex const c1(v1), c2(v2); \
        BOOST_TEST(c1.real() == c2.real()); \
        BOOST_TEST(c1.imag() == c2.imag()); \
    }

BOOST_AUTO_TEST_CASE(Empty)
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms);
        
    sycomore::epg::Regular model(species);
    BOOST_TEST(model.states_count() == 1);
    BOOST_TEST((model.state(0) == std::vector<sycomore::Complex>{0,0,1}));
    BOOST_TEST((model.states() == std::vector<sycomore::Complex>{0,0,1}));
    BOOST_TEST(model.echo() == 0.);
}

BOOST_AUTO_TEST_CASE(Pulse, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    
    BOOST_TEST(model.states_count() == 1);
    auto const state = model.state(0);
    
    TEST_COMPLEX_EQUAL(
        state[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(
        state[1], sycomore::Complex(0.2857626571584661, +0.6732146319308543));
    TEST_COMPLEX_EQUAL(state[2], 0.6819983600624985);

    std::vector<sycomore::Complex> const states{
        {0.2857626571584661, -0.6732146319308543},
        {0.2857626571584661, +0.6732146319308543},
        0.6819983600624985
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        TEST_COMPLEX_EQUAL(model.states()[i], states[i]);
    }
}

BOOST_AUTO_TEST_CASE(Gradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift();
    
    BOOST_TEST(model.states_count() == 2);
    
    auto const state_0 = model.state(0);
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6819983600624985);
    
    auto const state_1 = model.state(1);
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6819983600624985,
        {0.2857626571584661, -0.6732146319308543}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        TEST_COMPLEX_EQUAL(model.states()[i], states[i]);
    }
}

BOOST_AUTO_TEST_CASE(Relaxation, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift();
    model.relaxation(10*ms);
    
    BOOST_TEST(model.states_count() == 2);
    
    auto const state_0 = model.state(0);
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);
    
    auto const state_1 = model.state(1);
    
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2585687448743616, -0.6091497893403431));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6851625292479138,
        {0.2585687448743616, -0.6091497893403431}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        TEST_COMPLEX_EQUAL(model.states()[i], states[i]);
    }
}

BOOST_AUTO_TEST_CASE(Diffusion, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift();
    model.relaxation(10*ms);
    model.diffusion(10*ms, 2*mT/m);
    
    BOOST_TEST(model.states_count() == 2);
    
    auto const state_0 = model.state(0);
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);
    
    auto const state_1 = model.state(1);
    
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.25805111586158685, -0.6079303318059787));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6851625292479138,
        {0.25805111586158685, -0.6079303318059787}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        TEST_COMPLEX_EQUAL(model.states()[i], states[i]);
    }
}

BOOST_AUTO_TEST_CASE(TimeInterval, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 2*mT/m);
    
    BOOST_TEST(model.states_count() == 2);
    
    auto const state_0 = model.state(0);
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);
    
    auto const state_1 = model.state(1);
    
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2584947343504123, -0.6089754314724013));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6851625292479138,
        {0.2584947343504123, -0.6089754314724013}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        TEST_COMPLEX_EQUAL(model.states()[i], states[i]);
    }
}
