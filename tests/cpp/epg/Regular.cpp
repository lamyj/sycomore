#define BOOST_TEST_MODULE epg_Regular
#include <boost/test/unit_test.hpp>

#include <xtensor/xview.hpp>

#include "sycomore/epg/Regular.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

#define TEST_COMPLEX_EQUAL(v1, v2) \
    { \
        sycomore::Complex const c1(v1), c2(v2); \
        BOOST_TEST(c1.real() == c2.real()); \
        BOOST_TEST(c1.imag() == c2.imag()); \
    }

#define TEST_ORDER(o1, o2) \
    { \
        BOOST_TEST(o1.magnitude == o2.magnitude); \
        BOOST_TEST(o1.dimensions == o2.dimensions); \
    }


void test_model(
    sycomore::epg::Regular const & model,
    sycomore::TensorQ<1> const & expected_orders,
    sycomore::ArrayC const & expected_states)
{
    auto && orders = model.orders();
    BOOST_TEST(orders.shape() == expected_orders.shape());
    
    auto && states = model.states();
    BOOST_TEST(model.size() == expected_states.shape()[0]);
    BOOST_TEST(states.shape() == expected_states.shape());
    
    for(std::size_t i=0; i<model.size(); ++i)
    {
        auto && order = orders(i);
        auto && expected_order = expected_orders(i);
        
        TEST_ORDER(order, expected_order);
        
        auto && expected_state = xt::view(expected_states, i);
        
        for(std::size_t j=0; j<expected_state.size(); ++j)
        {
            TEST_COMPLEX_EQUAL(model.state(i)(j), expected_state(j));
            TEST_COMPLEX_EQUAL(model.state(order)(j), expected_state(j));
        }
    }
    BOOST_TEST(model.echo() == expected_states.at(0, 0));
}

BOOST_AUTO_TEST_CASE(Empty)
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms);
        
    sycomore::epg::Regular model(species);
    
    test_model(model, {0}, {{0,0,1}});
}

BOOST_AUTO_TEST_CASE(Pulse, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    
    test_model(
        model, 
        {0},
        {
            {
                {0.2857626571584661, -0.6732146319308543},
                {0.2857626571584661, +0.6732146319308543},
                0.6819983600624985}});
}

BOOST_AUTO_TEST_CASE(Gradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift();
    
    test_model(
        model,
        {0, 1},
        {
            {0, 0, 0.6819983600624985},
            {{0.2857626571584661, -0.6732146319308543}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(Relaxation, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift();
    model.relaxation(10*ms);
    
    test_model(
        model, 
        {0, 1},
        {
            {0, 0, 0.6851625292479138},
            {{0.2585687448743616, -0.6091497893403431}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(Diffusion, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species, {0,0,1}, 100, 20*mT/m*ms);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);
    model.relaxation(10*ms);
    model.diffusion(10*ms, 2*mT/m);
    
    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms, sycomore::gamma*20*mT/m*ms},
        {
            {0, 0, 0.6851625292479138},
            {{0.25805111586158685, -0.60793033180597855}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(OffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
    
    sycomore::epg::Regular model(species);
    model.delta_omega = 10*Hz;
    model.apply_pulse(47*deg, 23*deg);
    model.shift();

    model.off_resonance(10*ms);
    test_model(
        model, 
        {0, 1},
        {
            {0, 0, 0.6819983600624985},
            {{0.6268924782754024, -0.37667500256027975}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(TimeInterval, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species, {0,0,1}, 100, 20*mT/m*ms);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 2*mT/m);
    
    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms, sycomore::gamma*20*mT/m*ms},
        {
            {0, 0, 0.6851625292479138},
            {{0.2584947343504123, -0.6089754314724013}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(TimeIntervalFieldOffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Regular model(species, {0,0,1}, 100, 20*mT/m*ms);
    model.delta_omega = 10*Hz;
    model.apply_pulse(47*deg, 23*deg);

    model.apply_time_interval(10*ms, 2*mT/m);
    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms, sycomore::gamma*20*mT/m*ms},
        {
            {0, 0, 0.6851625292479138},
            {{0.56707341067384409, -0.34073208057155585}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(TimeIntervalSpeciesOffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms, 10*Hz);

    sycomore::epg::Regular model(species, {0,0,1}, 100, 20*mT/m*ms);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 2*mT/m);

    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms, sycomore::gamma*20*mT/m*ms},
        {
            {0, 0, 0.6851625292479138},
            {{0.56707341067384409, -0.34073208057155585}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(TimeIntervalBothOffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms, 10*Hz);

    sycomore::epg::Regular model(species, {0,0,1}, 100, 20*mT/m*ms);
    model.delta_omega = -model.species().delta_omega();
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval({10*ms, 2*mT/m});
    
    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms, sycomore::gamma*20*mT/m*ms},
        {
            {0, 0, 0.6851625292479138},
            {{0.2584947343504123, -0.6089754314724013}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(UnitGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms);
        
    sycomore::epg::Regular model(species, {0,0,1}, 100, 10*mT/m*ms);
    model.apply_pulse(47*deg, 23*deg);
    
    /* First time interval: 1*unit dephasing */
    model.apply_time_interval(10*ms, 1*mT/m);
    
    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms, sycomore::gamma*10*mT/m*ms},
        {
            {0, 0, 0.6851625292479138},
            {{0.2585687448743616, -0.609149789340343}, 0, 0}});
    
    /* Second time interval: 2*unit dephasing */
    model.apply_time_interval(10*ms, 2*mT/m);
    test_model(
        model, 
        {
            sycomore::gamma*0*mT/m*ms, sycomore::gamma*10*mT/m*ms,
            sycomore::gamma*20*mT/m*ms, sycomore::gamma*30*mT/m*ms},
        {
            {0, 0, 0.6882952144238884},
            {0, 0, 0},
            {0, 0, 0},
            {{0.2339626754969161, -0.5511815225838647}, 0, 0}});
    
    /* Third time interval: -3*unit dephasing. Higher states are all 0. */
    model.apply_time_interval(10*ms, -3*mT/m);
    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms},
        {
            {
                {0.2116981832134146, -0.49872966576391303}, 
                {0.2116981832134146, +0.49872966576391303}, 
                0.6913967288615507},
        });
    
    BOOST_CHECK_THROW(
        model.apply_time_interval(12*ms, 2*mT/m), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(BulkMotion, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms);
        
    sycomore::epg::Regular model(species, {0,0,1}, 100, 10*mT/m*ms);
    model.velocity = 40*cm/s;
    
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 1*mT/m);
    
    test_model(
        model, 
        {sycomore::gamma*0*mT/m*ms, sycomore::gamma*10*mT/m*ms},
        {
            {0, 0, 0.6851625292479138},
            {{-0.33529082747796918, -0.57052723220581303}, 0, 0}});
}

BOOST_AUTO_TEST_CASE(Elapsed)
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms);
        
    sycomore::epg::Regular model(species, {0,0,1}, 100, 10*mT/m*ms);
    BOOST_TEST(model.elapsed() == 0*s);
    
    model.apply_time_interval(10*ms);
    BOOST_TEST(model.elapsed() == 10*ms);
}
