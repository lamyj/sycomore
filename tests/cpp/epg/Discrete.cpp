#define BOOST_TEST_MODULE epg_Discrete
#include <boost/test/unit_test.hpp>

#include "sycomore/epg/Discrete.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(Empty)
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms);
        
    sycomore::epg::Discrete model(species);
    BOOST_TEST((model.orders() == std::vector<sycomore::Quantity>{0*rad/m}));
    BOOST_TEST((model.state(0) == std::vector<sycomore::Complex>{0,0,1}));
    BOOST_TEST(model.echo() == 0.);
}

BOOST_AUTO_TEST_CASE(Pulse, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);

    BOOST_TEST((model.orders() == std::vector<sycomore::Quantity>{0*rad/m}));
    auto const state = model.state(0);

    BOOST_TEST(state[0].real() == 0.2857626571584661);
    BOOST_TEST(state[0].imag() == -0.6732146319308543);

    BOOST_TEST(state[1].real() == 0.2857626571584661);
    BOOST_TEST(state[1].imag() == +0.6732146319308543);

    BOOST_TEST(state[2].real() == 0.6819983600624985);
    BOOST_TEST(state[2].imag() == 0.);

    std::vector<sycomore::Complex> const states{
        {0.2857626571584661, -0.6732146319308543},
        {0.2857626571584661, +0.6732146319308543},
        0.6819983600624985
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        BOOST_TEST(model.states()[i].real() == states[i].real());
        BOOST_TEST(model.states()[i].imag() == states[i].imag());
    }
}

BOOST_AUTO_TEST_CASE(PositiveGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);
    
    BOOST_TEST((
        model.orders() == std::vector<sycomore::Quantity>{0*rad/m, 5350*rad/m}));

    auto const state_0 = model.state(0);
    BOOST_TEST(state_0[0] == 0.);
    BOOST_TEST(state_0[1] == 0.);
    BOOST_TEST(state_0[2] == 0.6819983600624985);

    auto const state_1 = model.state(1);
    BOOST_TEST(state_1[0].real() == 0.2857626571584661);
    BOOST_TEST(state_1[0].imag() == -0.6732146319308543);
    BOOST_TEST(state_1[1] == 0.);
    BOOST_TEST(state_1[2] == 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6819983600624985,
        {0.2857626571584661, -0.6732146319308543}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        BOOST_TEST(model.states()[i].real() == states[i].real());
        BOOST_TEST(model.states()[i].imag() == states[i].imag());
    }
}

BOOST_AUTO_TEST_CASE(NegativeGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, -2*mT/m);
    
    BOOST_TEST((
        model.orders() == std::vector<sycomore::Quantity>{0*rad/m, 5350*rad/m}));
    
    auto const state_0 = model.state(0);
    BOOST_TEST(state_0[0] == 0.);
    BOOST_TEST(state_0[1] == 0.);
    BOOST_TEST(state_0[2] == 0.6819983600624985);

    auto const state_1 = model.state(1);
    BOOST_TEST(state_1[0] == 0.);
    BOOST_TEST(state_1[1].real() == 0.2857626571584661);
    BOOST_TEST(state_1[1].imag() == 0.6732146319308543);
    BOOST_TEST(state_1[2] == 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6819983600624985,
        0, {0.2857626571584661, 0.6732146319308543}, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        BOOST_TEST(model.states()[i].real() == states[i].real());
        BOOST_TEST(model.states()[i].imag() == states[i].imag());
    }
}

BOOST_AUTO_TEST_CASE(MultipleGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, -2*mT/m);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 1*mT/m);
    
    BOOST_TEST((
        model.orders() == std::vector<sycomore::Quantity>{
            0*rad/m, 2675*rad/m, 5350*rad/m, 8025*rad/m}));
    
    auto const state_0 = model.state(0);
    BOOST_TEST(state_0[0] == 0.);
    BOOST_TEST(state_0[1] == 0.);
    BOOST_TEST(state_0[2] == 0.4651217631279373);

    auto const state_1 = model.state(1);
    BOOST_TEST(state_1[0].real() == 0.19488966354917586);
    BOOST_TEST(state_1[0].imag() == -0.45913127494692113);
    BOOST_TEST(state_1[1].real() == 0.240326160353821);
    BOOST_TEST(state_1[1].imag() == 0.5661729534388877);
    BOOST_TEST(state_1[2] == 0.);
    
    auto const state_2 = model.state(2);
    BOOST_TEST(state_2[0] == 0.);
    BOOST_TEST(state_2[1] == 0.);
    BOOST_TEST(state_2[2] == -0.26743911843603135);
    
    auto const state_3 = model.state(3);
    BOOST_TEST(state_3[0].real() == -0.045436496804645087);
    BOOST_TEST(state_3[0].imag() == 0.10704167849196657);
    BOOST_TEST(state_3[1] == 0.);
    BOOST_TEST(state_3[2] == 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.4651217631279373,

        {0.19488966354917586, -0.45913127494692113},
        {0.240326160353821, 0.5661729534388877},
        0,

        0, 0, -0.26743911843603135,

        {-0.045436496804645087, 0.10704167849196657}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        BOOST_TEST(model.states()[i].real() == states[i].real());
        BOOST_TEST(model.states()[i].imag() == states[i].imag());
    }
}


BOOST_AUTO_TEST_CASE(Relaxation, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);
    model.relaxation(10*ms);

    BOOST_TEST((
        model.orders() == std::vector<sycomore::Quantity>{0*rad/m, 5350*rad/m}));

    auto const state_0 = model.state(0);
    BOOST_TEST(state_0[0] == 0.);
    BOOST_TEST(state_0[1] == 0.);
    BOOST_TEST(state_0[2] == 0.6851625292479138);

    auto const state_1 = model.state(1);

    BOOST_TEST(state_1[0].real() == 0.2585687448743616);
    BOOST_TEST(state_1[0].imag() == -0.6091497893403431);
    BOOST_TEST(state_1[1] == 0.);
    BOOST_TEST(state_1[2] == 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6851625292479138,
        {0.2585687448743616, -0.6091497893403431}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        BOOST_TEST(model.states()[i].real() == states[i].real());
        BOOST_TEST(model.states()[i].imag() == states[i].imag());
    }
}

BOOST_AUTO_TEST_CASE(Diffusion, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, 2*mT/m);
    model.relaxation(10*ms);
    model.diffusion(10*ms, 2*mT/m);

    BOOST_TEST((
        model.orders() == std::vector<sycomore::Quantity>{0*rad/m, 5350*rad/m}));

    auto const state_0 = model.state(0);
    BOOST_TEST(state_0[0] == 0.);
    BOOST_TEST(state_0[1] == 0.);
    BOOST_TEST(state_0[2] == 0.6851625292479138);

    auto const state_1 = model.state(1);

    BOOST_TEST(state_1[0].real() == 0.25805117100742553);
    BOOST_TEST(state_1[0].imag() == -0.6079304617214332);
    BOOST_TEST(state_1[1] == 0.);
    BOOST_TEST(state_1[2] == 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6851625292479138,
        {0.25805117100742553, -0.6079304617214332}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        BOOST_TEST(model.states()[i].real() == states[i].real());
        BOOST_TEST(model.states()[i].imag() == states[i].imag());
    }
}

BOOST_AUTO_TEST_CASE(TimeInterval, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 2*mT/m);

    BOOST_TEST((
        model.orders() == std::vector<sycomore::Quantity>{0*rad/m, 5350*rad/m}));

    auto const state_0 = model.state(0);
    BOOST_TEST(state_0[0] == 0.);
    BOOST_TEST(state_0[1] == 0.);
    BOOST_TEST(state_0[2] == 0.6851625292479138);

    auto const state_1 = model.state(1);

    BOOST_TEST(state_1[0].real() == 0.2584947343504123);
    BOOST_TEST(state_1[0].imag() == -0.6089754314724013);
    BOOST_TEST(state_1[1] == 0.);
    BOOST_TEST(state_1[2] == 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6851625292479138,
        {0.2584947343504123, -0.6089754314724013}, 0, 0
    };
    BOOST_TEST(model.states().size() == states.size());
    for(int i=0; i<model.states().size(); ++i)
    {
        BOOST_TEST(model.states()[i].real() == states[i].real());
        BOOST_TEST(model.states()[i].imag() == states[i].imag());
    }
}
