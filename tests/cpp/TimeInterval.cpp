#define BOOST_TEST_MODULE TimeInterval
#include <boost/test/unit_test.hpp>

#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

void test_quantity_array(
    sycomore::Vector3<sycomore::Quantity> const & left,
    sycomore::Vector3<sycomore::Quantity> const & right)
{
    BOOST_TEST(left.size() == right.size());
    for(unsigned int i=0; i<left.size(); ++i)
    {
        auto const & l = left[i];
        auto const & r = right[i];
        
        BOOST_TEST(l.dimensions == r.dimensions);
        BOOST_TEST(l.magnitude == r.magnitude);
    }
}

sycomore::Vector3<sycomore::Quantity> const amplitude {
    20*sycomore::units::mT/sycomore::units::m, 
    40*sycomore::units::mT/sycomore::units::m, 
    80*sycomore::units::mT/sycomore::units::m};
sycomore::Vector3<sycomore::Quantity> const area(1*sycomore::units::ms*amplitude);
sycomore::Vector3<sycomore::Quantity> const dephasing(sycomore::gamma*area);

BOOST_AUTO_TEST_CASE(DefaultConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval;
    
    BOOST_TEST(interval.duration() == 0*s);
    
    test_quantity_array(interval.gradient_amplitude(), {0*T/m, 0*T/m, 0*T/m});
    test_quantity_array(interval.gradient_area(), {0*T/m*s, 0*T/m*s, 0*T/m*s});
    test_quantity_array(
        interval.gradient_dephasing(), {0*rad/m, 0*rad/m, 0*rad/m});
}

BOOST_AUTO_TEST_CASE(DurationConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms);
    
    BOOST_TEST(interval.duration() == 1*ms);
    
    test_quantity_array(interval.gradient_amplitude(), {0*T/m, 0*T/m, 0*T/m});
    test_quantity_array(interval.gradient_area(), {0*T/m*s, 0*T/m*s, 0*T/m*s});
    test_quantity_array(
        interval.gradient_dephasing(), {0*rad/m, 0*rad/m, 0*rad/m});
}

BOOST_AUTO_TEST_CASE(DephasingScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, dephasing[0]);
    
    BOOST_TEST(interval.duration() == 1*ms);
    
    test_quantity_array(
        interval.gradient_amplitude(),
        {amplitude[0], amplitude[0], amplitude[0]});
    test_quantity_array(interval.gradient_area(), {area[0], area[0], area[0]});
    test_quantity_array(
        interval.gradient_dephasing(),
        {dephasing[0], dephasing[0], dephasing[0]});
}

BOOST_AUTO_TEST_CASE(DephasingVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, dephasing);
    
    BOOST_TEST(interval.duration() == 1*ms);
    
    test_quantity_array(interval.gradient_amplitude(), amplitude);
    test_quantity_array(interval.gradient_area(), area);
    test_quantity_array(interval.gradient_dephasing(), dephasing);
}

BOOST_AUTO_TEST_CASE(AmplitudeScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, amplitude[0]);
    
    BOOST_TEST(interval.duration() == 1*ms);
    
    test_quantity_array(
        interval.gradient_amplitude(),
        {amplitude[0], amplitude[0], amplitude[0]});
    test_quantity_array(interval.gradient_area(), {area[0], area[0], area[0]});
    test_quantity_array(
        interval.gradient_dephasing(),
        {dephasing[0], dephasing[0], dephasing[0]});
}

BOOST_AUTO_TEST_CASE(AmplitudeVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, amplitude);
    
    BOOST_TEST(interval.duration() == 1*ms);
    
    test_quantity_array(interval.gradient_amplitude(), amplitude);
    test_quantity_array(interval.gradient_area(), area);
    test_quantity_array(interval.gradient_dephasing(), dephasing);
}

BOOST_AUTO_TEST_CASE(AreaScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, area[0]);
    
    BOOST_TEST(interval.duration() == 1*ms);
    
    test_quantity_array(
        interval.gradient_amplitude(),
        {amplitude[0], amplitude[0], amplitude[0]});
    test_quantity_array(interval.gradient_area(), {area[0], area[0], area[0]});
    test_quantity_array(
        interval.gradient_dephasing(),
        {dephasing[0], dephasing[0], dephasing[0]});
}

BOOST_AUTO_TEST_CASE(AreaVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1*ms, area);
    
    BOOST_TEST(interval.duration() == 1*ms);
    
    test_quantity_array(interval.gradient_amplitude(), amplitude);
    test_quantity_array(interval.gradient_area(), area);
    test_quantity_array(interval.gradient_dephasing(), dephasing);
}

void test_gradient_accessor(sycomore::Vector3<sycomore::Quantity> const & data)
{
    using namespace sycomore::units;
    sycomore::TimeInterval interval(1._ms);
    
    interval.set_gradient(data[1]);
    test_quantity_array(
        interval.gradient_amplitude(), 
        {amplitude[1], amplitude[1], amplitude[1]});
    test_quantity_array(interval.gradient_area(), {area[1], area[1], area[1]});
    test_quantity_array(
        interval.gradient_dephasing(), 
        {dephasing[1], dephasing[1], dephasing[1]});
    
    interval.set_gradient(data);
    test_quantity_array(interval.gradient_amplitude(), amplitude);
    test_quantity_array(interval.gradient_area(), area);
    test_quantity_array(interval.gradient_dephasing(), dephasing);
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
    BOOST_TEST(interval_1.gradient_area()[0] == 100*mT/m*ms);
    BOOST_TEST(interval_1.gradient_amplitude()[0] == G_max);
    
    auto const interval_2 = sycomore::TimeInterval::shortest(
        sycomore::gamma*100*mT/m*ms, G_max);
    BOOST_TEST(
        interval_2.gradient_area()[0].magnitude == (100*mT/m*ms).magnitude);
    BOOST_TEST(interval_2.gradient_amplitude()[0] == G_max);
}

BOOST_AUTO_TEST_CASE(Shortest_3D, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    
    auto const G_max = 20*mT/m;
    
    auto const interval_1 = sycomore::TimeInterval::shortest(
        {100*mT/m*ms, 200*mT/m*ms, 400*mT/m*ms}, G_max);
    BOOST_TEST(interval_1.gradient_area()[0] == 100*mT/m*ms);
    BOOST_TEST(interval_1.gradient_area()[1] == 200*mT/m*ms);
    BOOST_TEST(interval_1.gradient_area()[2] == 400*mT/m*ms);
    
    BOOST_TEST(interval_1.gradient_amplitude()[0] < G_max);
    BOOST_TEST(interval_1.gradient_amplitude()[1] < G_max);
    BOOST_TEST(interval_1.gradient_amplitude()[2] == G_max);
    
    auto const interval_2 = sycomore::TimeInterval::shortest(
        sycomore::gamma*sycomore::Vector3Q{100*mT/m*ms,200*mT/m*ms,400*mT/m*ms},
        G_max);
    BOOST_TEST(
        interval_2.gradient_area()[0].magnitude == (100*mT/m*ms).magnitude);
    BOOST_TEST(
        interval_2.gradient_area()[1].magnitude == (200*mT/m*ms).magnitude);
    BOOST_TEST(
        interval_2.gradient_area()[2].magnitude == (400*mT/m*ms).magnitude);
    BOOST_TEST(interval_2.gradient_amplitude()[0] < G_max);
    BOOST_TEST(interval_2.gradient_amplitude()[1] < G_max);
    BOOST_TEST(interval_2.gradient_amplitude()[2] == G_max);
}
