#define BOOST_TEST_MODULE TimeInterval
#include <boost/test/unit_test.hpp>

#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

#include "utils.h"

sycomore::Vector3Q const amplitude {
    20*sycomore::units::mT/sycomore::units::m, 
    40*sycomore::units::mT/sycomore::units::m, 
    80*sycomore::units::mT/sycomore::units::m};
sycomore::Vector3Q const area(1*sycomore::units::ms*amplitude);
sycomore::Vector3Q const dephasing(sycomore::gamma*area);

BOOST_AUTO_TEST_CASE(DefaultConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval;
    
    CHECK_QUANTITY(interval.duration(), 0*s);
    
    CHECK_QUANTITY(
        interval.gradient_amplitude(), (sycomore::Vector3Q{0*T/m, 0*T/m, 0*T/m}));
    CHECK_QUANTITY(
        interval.gradient_area(), (sycomore::Vector3Q{0*T/m*s, 0*T/m*s, 0*T/m*s}));
    CHECK_QUANTITY(
        interval.gradient_dephasing(),
        (sycomore::Vector3Q{0*rad/m, 0*rad/m, 0*rad/m}));
}

BOOST_AUTO_TEST_CASE(DurationConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms);
    
    CHECK_QUANTITY(interval.duration(), 1*ms);
    
    CHECK_QUANTITY(
        interval.gradient_amplitude(), (sycomore::Vector3Q{0*T/m, 0*T/m, 0*T/m}));
    CHECK_QUANTITY(
        interval.gradient_area(), (sycomore::Vector3Q{0*T/m*s, 0*T/m*s, 0*T/m*s}));
    CHECK_QUANTITY(
        interval.gradient_dephasing(), (sycomore::Vector3Q{0*rad/m, 0*rad/m, 0*rad/m}));
}

BOOST_AUTO_TEST_CASE(DephasingScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, dephasing[0]);
    
    CHECK_QUANTITY(interval.duration(), 1*ms);
    
    CHECK_QUANTITY(
        interval.gradient_amplitude(),
        (sycomore::Vector3Q{amplitude[0], amplitude[0], amplitude[0]}));
    CHECK_QUANTITY(
        interval.gradient_area(),
        (sycomore::Vector3Q{area[0], area[0], area[0]}));
    CHECK_QUANTITY(
        interval.gradient_dephasing(),
        (sycomore::Vector3Q{dephasing[0], dephasing[0], dephasing[0]}));
}

BOOST_AUTO_TEST_CASE(DephasingVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, dephasing);
    
    CHECK_QUANTITY(interval.duration(), 1*ms);
    
    CHECK_QUANTITY(interval.gradient_amplitude(), amplitude);
    CHECK_QUANTITY(interval.gradient_area(), area);
    CHECK_QUANTITY(interval.gradient_dephasing(), dephasing);
}

BOOST_AUTO_TEST_CASE(AmplitudeScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, amplitude[0]);
    
    CHECK_QUANTITY(interval.duration(), 1*ms);
    
    CHECK_QUANTITY(
        interval.gradient_amplitude(),
        (sycomore::Vector3Q{amplitude[0], amplitude[0], amplitude[0]}));
    CHECK_QUANTITY(
        interval.gradient_area(),
        (sycomore::Vector3Q{area[0], area[0], area[0]}));
    CHECK_QUANTITY(
        interval.gradient_dephasing(),
        (sycomore::Vector3Q{dephasing[0], dephasing[0], dephasing[0]}));
}

BOOST_AUTO_TEST_CASE(AmplitudeVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, amplitude);
    
    CHECK_QUANTITY(interval.duration(), 1*ms);
    
    CHECK_QUANTITY(interval.gradient_amplitude(), amplitude);
    CHECK_QUANTITY(interval.gradient_area(), area);
    CHECK_QUANTITY(interval.gradient_dephasing(), dephasing);
}

BOOST_AUTO_TEST_CASE(AreaScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, area[0]);
    
    CHECK_QUANTITY(interval.duration(), 1*ms);
    
    CHECK_QUANTITY(
        interval.gradient_amplitude(),
        (sycomore::Vector3Q{amplitude[0], amplitude[0], amplitude[0]}));
    CHECK_QUANTITY(
        interval.gradient_area(),
        (sycomore::Vector3Q{area[0], area[0], area[0]}));
    CHECK_QUANTITY(
        interval.gradient_dephasing(),
        (sycomore::Vector3Q{dephasing[0], dephasing[0], dephasing[0]}));
}

BOOST_AUTO_TEST_CASE(AreaVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, area);
    
    CHECK_QUANTITY(interval.duration(), 1*ms);
    
    CHECK_QUANTITY(interval.gradient_amplitude(), amplitude);
    CHECK_QUANTITY(interval.gradient_area(), area);
    CHECK_QUANTITY(interval.gradient_dephasing(), dephasing);
}

void test_gradient_accessor(sycomore::Vector3Q const & data)
{
    using namespace sycomore::units;
    sycomore::TimeInterval interval(1._ms);
    
    interval.set_gradient(data[1]);
    CHECK_QUANTITY(
        interval.gradient_amplitude(), 
        (sycomore::Vector3Q{amplitude[1], amplitude[1], amplitude[1]}));
    CHECK_QUANTITY(
        interval.gradient_area(),
        (sycomore::Vector3Q{area[1], area[1], area[1]}));
    CHECK_QUANTITY(
        interval.gradient_dephasing(), 
        (sycomore::Vector3Q{dephasing[1], dephasing[1], dephasing[1]}));
    
    interval.set_gradient(data);
    CHECK_QUANTITY(interval.gradient_amplitude(), amplitude);
    CHECK_QUANTITY(interval.gradient_area(), area);
    CHECK_QUANTITY(interval.gradient_dephasing(), dephasing);
}

BOOST_AUTO_TEST_CASE(GradientAccessors, *boost::unit_test::tolerance(1e-9))
{
    test_gradient_accessor(amplitude);
    test_gradient_accessor(area);
    test_gradient_accessor(dephasing);
}

BOOST_AUTO_TEST_CASE(Comparison)
{
    using namespace sycomore::units;
    
    sycomore::TimeInterval const interval_1(1._ms, 2*T/m);
    sycomore::TimeInterval const interval_2(1._ms, 2*T/m);
    sycomore::TimeInterval const interval_3(1._ms, {2*T/m, 2*T/m, 2*T/m});
    
    BOOST_CHECK(interval_1 == interval_2);
    BOOST_CHECK(interval_1 == interval_3);
    BOOST_CHECK(!(interval_1 != interval_2));
    BOOST_CHECK(!(interval_1 != interval_3));
    
    sycomore::TimeInterval const interval_4(4._ms, 2*T/m);
    BOOST_CHECK(!(interval_1 == interval_4));
    BOOST_CHECK(interval_1 != interval_4);
    
    sycomore::TimeInterval const interval_5(1._ms, 4*T/m);
    BOOST_CHECK(!(interval_1 == interval_5));
    BOOST_CHECK(interval_1 != interval_5);
}

BOOST_AUTO_TEST_CASE(Shortest_1D, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    
    auto const G_max = 20*mT/m;
    
    auto const interval_1 = sycomore::TimeInterval::shortest(
        100*mT/m*ms, G_max);
    CHECK_QUANTITY(interval_1.gradient_area()[0], 100*mT/m*ms);
    CHECK_QUANTITY(interval_1.gradient_amplitude()[0], G_max);
    
    auto const interval_2 = sycomore::TimeInterval::shortest(
        sycomore::gamma*100*mT/m*ms, G_max);
    CHECK_QUANTITY(interval_2.gradient_area()[0], 100*mT/m*ms);
    CHECK_QUANTITY(interval_2.gradient_amplitude()[0], G_max);
}

BOOST_AUTO_TEST_CASE(Shortest_3D, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    
    auto const G_max = 20*mT/m;
    
    auto const interval_1 = sycomore::TimeInterval::shortest(
        {100*mT/m*ms, 200*mT/m*ms, 400*mT/m*ms}, G_max);
    CHECK_QUANTITY(interval_1.gradient_area()[0], 100*mT/m*ms);
    CHECK_QUANTITY(interval_1.gradient_area()[1], 200*mT/m*ms);
    CHECK_QUANTITY(interval_1.gradient_area()[2], 400*mT/m*ms);
    
    BOOST_TEST((interval_1.gradient_amplitude()[0] < G_max));
    BOOST_TEST((interval_1.gradient_amplitude()[1] < G_max));
    CHECK_QUANTITY(interval_1.gradient_amplitude()[2], G_max);
    
    auto const interval_2 = sycomore::TimeInterval::shortest(
        sycomore::gamma*sycomore::Vector3Q{100*mT/m*ms,200*mT/m*ms,400*mT/m*ms},
        G_max);
    CHECK_QUANTITY(interval_2.gradient_area()[0], 100*mT/m*ms);
    CHECK_QUANTITY(interval_2.gradient_area()[1], 200*mT/m*ms);
    CHECK_QUANTITY(interval_2.gradient_area()[2], 400*mT/m*ms);
    BOOST_TEST((interval_2.gradient_amplitude()[0] < G_max));
    BOOST_TEST((interval_2.gradient_amplitude()[1] < G_max));
    CHECK_QUANTITY(interval_2.gradient_amplitude()[2], G_max);
}
