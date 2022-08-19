#define BOOST_TEST_MODULE epg_Discrete
#include <boost/test/unit_test.hpp>

#include "sycomore/epg/Discrete.h"
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
    sycomore::epg::Discrete const & model,
    std::vector<sycomore::epg::Discrete::Order> const & expected_orders,
    std::vector<sycomore::epg::Discrete::State> const & expected_states)
{
    auto && orders = model.orders();
    BOOST_TEST(orders.size() == expected_orders.size());

    auto && states = model.states();
    BOOST_TEST(model.size() == expected_states.size());
    BOOST_TEST(states.size() == 3*expected_states.size());

    for(std::size_t i=0; i<model.size(); ++i)
    {
        auto && expected_order = expected_orders[i];
        auto && order = orders[i];
        TEST_ORDER(order, expected_order);

        auto && expected_state = expected_states[i];
        {
            auto && state = model.state(i);
            BOOST_TEST(state.size() == expected_state.size());
            for(std::size_t j=0; j<state.size(); ++j)
            {
                TEST_COMPLEX_EQUAL(state[j], expected_state[j]);
            }
        }
        {
            auto && state = model.state(order);
            BOOST_TEST(state.size() == expected_state.size());
            for(std::size_t j=0; j<state.size(); ++j)
            {
                TEST_COMPLEX_EQUAL(state[j], expected_state[j]);
            }
        }
        {
            sycomore::epg::Discrete::State const state{
                states[3*i+0], states[3*i+1], states[3*i+2]};
            BOOST_TEST(state.size() == expected_state.size());
            for(std::size_t j=0; j<state.size(); ++j)
            {
                TEST_COMPLEX_EQUAL(state[j], expected_state[j]);
            }
        }
    }
}

sycomore::Species const species(
    1000*sycomore::units::ms, 100*sycomore::units::ms,
    3*sycomore::units::um*sycomore::units::um/sycomore::units::ms);

BOOST_AUTO_TEST_CASE(Empty, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
        
    sycomore::epg::Discrete model(species);

    std::vector<sycomore::epg::Discrete::Order> const orders{0*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{{0, 0, 1}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0.);
}

BOOST_AUTO_TEST_CASE(Pulse, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);

    std::vector<sycomore::epg::Discrete::Order> const orders{0*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {
            {0.2857626571584661, -0.6732146319308543},
            {0.2857626571584661, +0.6732146319308543},
            0.6819983600624985}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(
        model.echo(),
        sycomore::Complex(0.2857626571584661, -0.6732146319308543));
}

BOOST_AUTO_TEST_CASE(PositiveGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6819983600624985},
        {{0.2857626571584661, -0.6732146319308543}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(NegativeGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, -2*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6819983600624985},
        {0, {0.2857626571584661, 0.6732146319308543}, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(MultipleGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, -2*mT/m);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 1*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 2675*rad/m, 5350*rad/m, 8025*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.4651217631279373},
        {
            {0.19488966354917586, -0.45913127494692113},
            {0.240326160353821, 0.5661729534388877},
            0},
        {0, 0, -0.26743911843603135},
        {{-0.045436496804645087, 0.10704167849196657}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(Relaxation, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);
    model.relaxation(10*ms);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6851625292479138},
        {{0.2585687448743616, -0.6091497893403431}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(Diffusion, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);
    model.relaxation(10*ms);
    model.diffusion(10*ms, 2*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6851625292479138},
        {{0.25805117100742553, -0.6079304617214332}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(OffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    
    sycomore::epg::Discrete model(species);
    model.delta_omega = 10*Hz;
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);

    model.off_resonance(10*ms);
    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6819983600624985},
        {{0.6268924782754024, -0.37667500256027975}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(TimeInterval, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 2*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6851625292479138},
        {{0.2584947343504123, -0.6089754314724013}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(TimeIntervalFieldOffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.delta_omega = 10*Hz;
    model.apply_pulse(47*deg, 23*deg);

    model.apply_time_interval({10*ms, 2*mT/m});
    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6851625292479138},
        {{0.56707341067384409, -0.34073208057155585}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(TimeIntervalSpeciesOffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(
        {species.get_R1(), species.get_R2(), species.get_D(), 0*Hz, 10*Hz});
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 2*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6851625292479138},
        {{0.56707341067384409, -0.34073208057155585}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(TimeIntervalBothOffResonance, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(
        {species.get_R1(), species.get_R2(), species.get_D(), 0*Hz, 10*Hz});
    model.delta_omega = -model.get_species().get_delta_omega();
    model.apply_pulse(47*deg, 23*deg);

    model.apply_time_interval({10*ms, 2*mT/m, });
    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {0, 0, 0.6851625292479138},
        {{0.2584947343504123, -0.6089754314724013}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(model.echo(), 0);
}

BOOST_AUTO_TEST_CASE(Refocalization, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::epg::Discrete model(species);
    model.apply_pulse(90*deg, 30*deg);
    model.apply_time_interval(10*ms, 2*mT/m);
    model.apply_pulse(120*deg, 0*deg);
    model.apply_time_interval(10*ms, 2*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 5350*rad/m, 10700*rad/m};
    std::vector<sycomore::epg::Discrete::State> const states{
        {
            {0.30684831950624042, 0.53147687960193668},
            {0.30684831950624042, -0.53147687960193668},
            0.0050245860296255166},
        {
            {0, -0.0077948398021822725},
            0,
            {-0.33555338970217136, -0.19373183987203996}},
        {{0.10210725404661349, -0.17685495183007738}, 0, 0}};
    test_model(model, orders, states);
    TEST_COMPLEX_EQUAL(
        model.echo(),
        sycomore::Complex(0.30684831950624042, 0.53147687960193668));
}

BOOST_AUTO_TEST_CASE(BulkMotion, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms);
        
    sycomore::epg::Discrete model(species);
    model.velocity = 40*cm/s;
    
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 1*mT/m);

    std::vector<sycomore::epg::Discrete::Order> const orders{
        0*rad/m, 2675*rad/m};    
    std::vector<sycomore::epg::Discrete::State> const states = {
        {0, 0, 0.6851625292479138},
        {{-0.33529079864474892, -0.57052724915068997}, 0, 0}};
    
    test_model(
        model, orders, states);
}
