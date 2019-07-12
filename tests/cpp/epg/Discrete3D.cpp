#define BOOST_TEST_MODULE epg_Discrete
#include <boost/test/unit_test.hpp>

#include "sycomore/epg/Discrete3D.h"
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

    sycomore::epg::Discrete3D model(species);
    BOOST_TEST(model.states().size() == 1);
    BOOST_TEST((model.state({0*rad/m,0*rad/m,0*rad/m}) == std::vector<sycomore::Complex>{0,0,1}));
    BOOST_TEST(model.echo() == 0.);
}

BOOST_AUTO_TEST_CASE(Pulse, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);

    BOOST_TEST(model.states().size() == 1);
    auto const state = model.state({0*rad/m,0*rad/m,0*rad/m});

    TEST_COMPLEX_EQUAL(
        state[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(
        state[1], sycomore::Complex(0.2857626571584661, +0.6732146319308543));
    TEST_COMPLEX_EQUAL(state[2], 0.6819983600624985);
}

BOOST_AUTO_TEST_CASE(PositiveGradientX, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {2*mT/m, 0*mT/m, 0*mT/m});

    BOOST_TEST(model.states().size() == 2);
    auto const state_0 = model.state({0*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6819983600624985);

    auto const state_1 = model.state({5350*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(state_1[1], sycomore::Complex(0.));
    TEST_COMPLEX_EQUAL(state_1[2], sycomore::Complex(0.));
}

BOOST_AUTO_TEST_CASE(PositiveGradientY, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {0*mT/m, 2*mT/m, 0*mT/m});

    BOOST_TEST(model.states().size() == 2);
    auto const state_0 = model.state({0*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6819983600624985);

    auto const state_1 = model.state({0*rad/m,5350*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(state_1[1], sycomore::Complex(0.));
    TEST_COMPLEX_EQUAL(state_1[2], sycomore::Complex(0.));
}

BOOST_AUTO_TEST_CASE(PositiveGradientZ, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {0*mT/m, 0*mT/m, 2*mT/m});

    BOOST_TEST(model.states().size() == 2);
    auto const state_0 = model.state({0*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6819983600624985);

    auto const state_1 = model.state({0*rad/m,0*rad/m,5350*rad/m});
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(state_1[1], sycomore::Complex(0.));
    TEST_COMPLEX_EQUAL(state_1[2], sycomore::Complex(0.));
}

BOOST_AUTO_TEST_CASE(NegativeGradientX, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {-2*mT/m, 0*mT/m, 0*mT/m});

    BOOST_TEST(model.states().size() == 2);
    auto const state_0 = model.state({0*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6819983600624985);

    auto const state_1 = model.state({5350*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(state_1[0], sycomore::Complex(0.));
    TEST_COMPLEX_EQUAL(
        state_1[1], sycomore::Complex(0.2857626571584661, 0.6732146319308543));
    TEST_COMPLEX_EQUAL(state_1[2], sycomore::Complex(0.));
}

BOOST_AUTO_TEST_CASE(NegativeGradientY, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {0*mT/m, -2*mT/m, 0*mT/m});

    BOOST_TEST(model.states().size() == 2);
    auto const state_0 = model.state({0*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6819983600624985);

    auto const state_1 = model.state({0*rad/m,-5350*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(state_1[1], sycomore::Complex(0.));
    TEST_COMPLEX_EQUAL(state_1[2], sycomore::Complex(0.));
}

BOOST_AUTO_TEST_CASE(NegativeGradientZ, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {0*mT/m, 0*mT/m, -2*mT/m});

    BOOST_TEST(model.states().size() == 2);
    auto const state_0 = model.state({0*rad/m,0*rad/m,0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6819983600624985);

    auto const state_1 = model.state({0*rad/m,0*rad/m,-5350*rad/m});
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2857626571584661, -0.6732146319308543));
    TEST_COMPLEX_EQUAL(state_1[1], sycomore::Complex(0.));
    TEST_COMPLEX_EQUAL(state_1[2], sycomore::Complex(0.));
}

BOOST_AUTO_TEST_CASE(MultipleGradient, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {-2*mT/m, 2*mT/m, -2*mT/m});
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {1*mT/m, -3*mT/m, 3*mT/m});

    BOOST_TEST(model.states().size() == 5);

    {
        auto const state = model.state({0*rad/m, 0*rad/m, 0*rad/m});
        TEST_COMPLEX_EQUAL(state[0], 0.);
        TEST_COMPLEX_EQUAL(state[1], 0.);
        TEST_COMPLEX_EQUAL(state[2], 0.4651217631279373);
    }

    {
        auto const state = model.state({2675*rad/m, 2676*rad/m, -2676*rad/m});
        TEST_COMPLEX_EQUAL(state[0], 0.);
        TEST_COMPLEX_EQUAL(
            state[1], sycomore::Complex(0.240326160353821, 0.5661729534388877));
        TEST_COMPLEX_EQUAL(state[2], 0.);
    }

    {
        auto const state = model.state({2675*rad/m, -8026*rad/m, 8026*rad/m});
        TEST_COMPLEX_EQUAL(
            state[0], sycomore::Complex(0.19488966354917586, -0.45913127494692113));
        TEST_COMPLEX_EQUAL(state[1], 0);
        TEST_COMPLEX_EQUAL(state[2], 0.);
    }

    {
        auto const state = model.state({5350*rad/m, -5350*rad/m, 5350*rad/m});
        TEST_COMPLEX_EQUAL(state[0], 0.);
        TEST_COMPLEX_EQUAL(state[1], 0.);
        TEST_COMPLEX_EQUAL(state[2], -0.26743911843603135);
    }

    {
        auto const state = model.state({8025*rad/m, -13376*rad/m, 13376*rad/m});
        TEST_COMPLEX_EQUAL(
            state[0], sycomore::Complex(-0.045436496804645087, 0.10704167849196657));
        TEST_COMPLEX_EQUAL(state[1], 0.);
        TEST_COMPLEX_EQUAL(state[2], 0.);
    }
}


BOOST_AUTO_TEST_CASE(Relaxation, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {2*mT/m, 0*mT/m, 0*mT/m});
    model.relaxation(10*ms);

    BOOST_TEST(model.states().size() == 2);

    auto const state_0 = model.state({0*rad/m, 0*rad/m, 0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);

    auto const state_1 = model.state({5350*rad/m, 0*rad/m, 0*rad/m});
    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2585687448743616, -0.6091497893403431));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);
}

BOOST_AUTO_TEST_CASE(DiffusionX, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {2*mT/m, 0*mT/m, 0*mT/m});
    model.relaxation(10*ms);
    model.diffusion(10*ms, {2*mT/m, 0*mT/m, 0*mT/m});

    BOOST_TEST(model.states().size() == 2);

    auto const state_0 = model.state({0*rad/m, 0*rad/m, 0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);

    auto const state_1 = model.state({5350*rad/m, 0*rad/m, 0*rad/m});

    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.25805117100742553, -0.6079304617214332));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);
}

BOOST_AUTO_TEST_CASE(DiffusionZ, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.shift(10*ms, {0*mT/m, 0*mT/m, 2*mT/m});
    model.relaxation(10*ms);
    model.diffusion(10*ms, {0*mT/m, 0*mT/m, 2*mT/m});

    BOOST_TEST(model.states().size() == 2);

    auto const state_0 = model.state({0*rad/m, 0*rad/m, 0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);

    auto const state_1 = model.state({0*rad/m, 0*rad/m, 5350*rad/m});

    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.25805117100742553, -0.6079304617214332));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);
}

BOOST_AUTO_TEST_CASE(TimeIntervalX, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, {2*mT/m, 0*mT/m, 0*mT/m});

    BOOST_TEST(model.states().size() == 2);

    auto const state_0 = model.state({0*rad/m, 0*rad/m, 0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);

    auto const state_1 = model.state({5350*rad/m, 0*rad/m, 0*rad/m});

    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2584947343504123, -0.6089754314724013));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);
}

BOOST_AUTO_TEST_CASE(TimeIntervalZ, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Species const species(1000*ms, 100*ms, 3*um*um/ms);

    sycomore::epg::Discrete3D model(species);
    model.apply_pulse(47*deg, 23*deg);
    model.apply_time_interval(10*ms, {0*mT/m, 0*mT/m, 2*mT/m});

    BOOST_TEST(model.states().size() == 2);

    auto const state_0 = model.state({0*rad/m, 0*rad/m, 0*rad/m});
    TEST_COMPLEX_EQUAL(state_0[0], 0.);
    TEST_COMPLEX_EQUAL(state_0[1], 0.);
    TEST_COMPLEX_EQUAL(state_0[2], 0.6851625292479138);

    auto const state_1 = model.state({0*rad/m, 0*rad/m, 5350*rad/m});

    TEST_COMPLEX_EQUAL(
        state_1[0], sycomore::Complex(0.2584947343504123, -0.6089754314724013));
    TEST_COMPLEX_EQUAL(state_1[1], 0.);
    TEST_COMPLEX_EQUAL(state_1[2], 0.);
}
