#define BOOST_TEST_MODULE TimeInterval
#include <boost/test/unit_test.hpp>

#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

void test_quantity_array(
    sycomore::Array<sycomore::Quantity> const & left,
    sycomore::Array<sycomore::Quantity> const & right)
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

BOOST_AUTO_TEST_CASE(DefaultConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval;
    
    BOOST_TEST(interval.get_duration() == 0*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {0*rad/m,0*rad/m,0*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {0*T/m,0*T/m,0*T/m});
    test_quantity_array(
        interval.get_gradient_area(), {0*T/m*s,0*T/m*s,0*T/m*s});
    test_quantity_array(
        interval.get_gradient_dephasing(), {0*rad/m,0*rad/m,0*rad/m});
}

BOOST_AUTO_TEST_CASE(DurationConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms);
    
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {0*rad/m,0*rad/m,0*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {0*T/m,0*T/m,0*T/m});
    test_quantity_array(
        interval.get_gradient_area(), {0*T/m*s,0*T/m*s,0*T/m*s});
    test_quantity_array(
        interval.get_gradient_dephasing(), {0*rad/m,0*rad/m,0*rad/m});
}

BOOST_AUTO_TEST_CASE(DephasingScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, 2.*rad/dm);
    
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {20*rad/m,20*rad/m,20*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {
            74.76015355016015*uT/m,
            74.76015355016015*uT/m,
            74.76015355016015*uT/m});
    test_quantity_array(
        interval.get_gradient_area(), {
            74.76015355016015*uT/m*ms,
            74.76015355016015*uT/m*ms,
            74.76015355016015*uT/m*ms});
    test_quantity_array(
        interval.get_gradient_dephasing(), {20*rad/m,20*rad/m,20*rad/m});
}

BOOST_AUTO_TEST_CASE(DephasingVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, {2*rad/dm,4*rad/m,8*rad/dam});
    
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {20*rad/m,4*rad/m,0.8*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {
            74.76015355016015*uT/m,
            14.95203071003203*uT/m,
            2.9904061420064063*uT/m});
    test_quantity_array(
        interval.get_gradient_area(), {
            74.76015355016015*uT/m*ms,
            14.95203071003203*uT/m*ms,
            2.9904061420064063*uT/m*ms});
    test_quantity_array(
        interval.get_gradient_dephasing(), {20*rad/m,4*rad/m,0.8*rad/m});
}

BOOST_AUTO_TEST_CASE(AmplitudeScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, 2.*mT/dm);
    
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {20*mT/m, 20*mT/m, 20*mT/m});
    test_quantity_array(
        interval.get_gradient_area(), {20*mT/m*ms, 20*mT/m*ms, 20*mT/m*ms});
    test_quantity_array(
        interval.get_gradient_dephasing(), {
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m});
}

BOOST_AUTO_TEST_CASE(AmplitudeVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, {2.*mT/dm, 4.*mT/m, 8.*mT/dam});
    
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {
            5350.4437993378515*rad/m,
            1070.0887598675702*rad/m,
            214.01775197351404*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {20*mT/m, 4*mT/m, 0.8*mT/m});
    test_quantity_array(
        interval.get_gradient_area(), {20*mT/m*ms, 4*mT/m*ms, 0.8*mT/m*ms});
    test_quantity_array(
        interval.get_gradient_dephasing(), {
            5350.4437993378515*rad/m,
            1070.0887598675702*rad/m,
            214.01775197351404*rad/m});
}

BOOST_AUTO_TEST_CASE(AreaScalarConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, 2.*mT/dm*ms);
    
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {20*mT/m, 20*mT/m, 20*mT/m});
    test_quantity_array(
        interval.get_gradient_area(), {20*mT/m*ms, 20*mT/m*ms, 20*mT/m*ms});
    test_quantity_array(
        interval.get_gradient_dephasing(), {
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m,
            5350.4437993378515*rad/m});
}

BOOST_AUTO_TEST_CASE(AreaVectorConstructor, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(
        1._ms, {2.*mT/dm*ms, 4.*mT/m*ms, 8.*mT/dam*ms});
    
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    
    test_quantity_array(
        interval.get_gradient_moment(), {
            5350.4437993378515*rad/m,
            1070.0887598675702*rad/m,
            214.01775197351404*rad/m});
    test_quantity_array(
        interval.get_gradient_amplitude(), {20*mT/m, 4*mT/m, 0.8*mT/m});
    test_quantity_array(
        interval.get_gradient_area(), {20*mT/m*ms, 4*mT/m*ms, 0.8*mT/m*ms});
    test_quantity_array(
        interval.get_gradient_dephasing(), {
            5350.4437993378515*rad/m,
            1070.0887598675702*rad/m,
            214.01775197351404*rad/m});
}

using Getter = sycomore::Quantity (sycomore::TimeInterval::*)(sycomore::Quantity const &);
using ScalarSetter = void (sycomore::TimeInterval::*)(sycomore::Quantity const &);
using VectorSetter = void (sycomore::TimeInterval::*)(sycomore::Array<sycomore::Quantity> const &);

sycomore::Array<sycomore::Quantity> const amplitude {
    20*sycomore::units::mT/sycomore::units::m, 
    40*sycomore::units::mT/sycomore::units::m, 
    80*sycomore::units::mT/sycomore::units::m};
sycomore::Array<sycomore::Quantity> const area {
    20*sycomore::units::mT/sycomore::units::m*sycomore::units::ms, 
    40*sycomore::units::mT/sycomore::units::m*sycomore::units::ms, 
    80*sycomore::units::mT/sycomore::units::m*sycomore::units::ms};
sycomore::Array<sycomore::Quantity> const dephasing{
    5350.4437993378515*sycomore::units::rad/sycomore::units::m,
    10700.887598675701*sycomore::units::rad/sycomore::units::m,
    21401.775197351402*sycomore::units::rad/sycomore::units::m};
auto const moment = dephasing;

template<typename TGetter>
void test_gradient_accessor(
    TGetter getter, ScalarSetter scalar_setter, VectorSetter vector_setter,
    sycomore::Array<sycomore::Quantity> const & data)
{
    using namespace sycomore::units;
    sycomore::TimeInterval interval(1._ms);
    
    (interval.*scalar_setter)(data[0]);
    test_quantity_array(
        interval.get_gradient_amplitude(), 
        {amplitude[0], amplitude[0], amplitude[0]});
    test_quantity_array(
        interval.get_gradient_area(), 
        {area[0], area[0], area[0]});
    test_quantity_array(
        interval.get_gradient_dephasing(), 
        {dephasing[0], dephasing[0], dephasing[0]});
    test_quantity_array(
        interval.get_gradient_moment(), 
        {dephasing[0], dephasing[0], dephasing[0]});
    
    (interval.*vector_setter)(data);
    test_quantity_array(interval.get_gradient_amplitude(), amplitude);
    test_quantity_array(interval.get_gradient_area(), area);
    test_quantity_array(interval.get_gradient_dephasing(), dephasing);
    test_quantity_array(interval.get_gradient_moment(), moment);
    
    interval.set_gradient(data[1]);
    test_quantity_array(
        interval.get_gradient_amplitude(), 
        {amplitude[1], amplitude[1], amplitude[1]});
    test_quantity_array(
        interval.get_gradient_area(), 
        {area[1], area[1], area[1]});
    test_quantity_array(
        interval.get_gradient_dephasing(), 
        {dephasing[1], dephasing[1], dephasing[1]});
    test_quantity_array(
        interval.get_gradient_moment(), 
        {dephasing[1], dephasing[1], dephasing[1]});
    
    interval.set_gradient(data);
    test_quantity_array(interval.get_gradient_amplitude(), amplitude);
    test_quantity_array(interval.get_gradient_area(), area);
    test_quantity_array(interval.get_gradient_dephasing(), dephasing);
    test_quantity_array(interval.get_gradient_moment(), moment);
}

BOOST_AUTO_TEST_CASE(GradientAccessors, *boost::unit_test::tolerance(1e-9))
{
    test_gradient_accessor(
        &sycomore::TimeInterval::get_gradient_amplitude,
        &sycomore::TimeInterval::set_gradient_amplitude,
        &sycomore::TimeInterval::set_gradient_amplitude,
        amplitude);
    test_gradient_accessor(
        &sycomore::TimeInterval::get_gradient_area,
        &sycomore::TimeInterval::set_gradient_area,
        &sycomore::TimeInterval::set_gradient_area,
        area);
    test_gradient_accessor(
        &sycomore::TimeInterval::get_gradient_dephasing,
        &sycomore::TimeInterval::set_gradient_dephasing,
        &sycomore::TimeInterval::set_gradient_dephasing,
        dephasing);
    test_gradient_accessor(
        &sycomore::TimeInterval::get_gradient_moment,
        &sycomore::TimeInterval::set_gradient_moment,
        &sycomore::TimeInterval::set_gradient_moment,
        moment);
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
