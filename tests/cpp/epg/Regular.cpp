#define BOOST_TEST_MODULE epg_Regular
#include <boost/test/unit_test.hpp>

#include "sycomore/epg/Regular.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

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
    
    BOOST_TEST(state[0].real() == 0.2857626571584661);
    BOOST_TEST(state[0].imag() == -0.6732146319308543);
    
    BOOST_TEST(state[1].real() == 0.2857626571584661);
    BOOST_TEST(state[1].imag() == +0.6732146319308543);
    
    BOOST_TEST(state[2].real() == 0.6819983600624985);
    BOOST_TEST(state[2].imag() == 0);

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

BOOST_AUTO_TEST_CASE(Gradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_gradient();
    
    BOOST_TEST(model.states_count() == 2);
    
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

BOOST_AUTO_TEST_CASE(Relaxation, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_gradient();
    model.apply_relaxation(10*ms);
    
    BOOST_TEST(model.states_count() == 2);
    
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
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_gradient();
    model.apply_relaxation(10*ms);
    model.apply_diffusion(10*ms, 2*mT/m);
    
    BOOST_TEST(model.states_count() == 2);
    
    auto const state_0 = model.state(0);
    BOOST_TEST(state_0[0] == 0.);
    BOOST_TEST(state_0[1] == 0.);
    BOOST_TEST(state_0[2] == 0.6851625292479138);
    
    auto const state_1 = model.state(1);
    
    BOOST_TEST(state_1[0].real() == 0.25805111586158685);
    BOOST_TEST(state_1[0].imag() == -0.6079303318059787);
    BOOST_TEST(state_1[1] == 0.);
    BOOST_TEST(state_1[2] == 0.);

    std::vector<sycomore::Complex> const states{
        0, 0, 0.6851625292479138,
        {0.25805111586158685, -0.6079303318059787}, 0, 0
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
        
    sycomore::epg::Regular model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, 2*mT/m);
    
    BOOST_TEST(model.states_count() == 2);
    
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
