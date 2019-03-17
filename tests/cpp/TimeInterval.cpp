#define BOOST_TEST_MODULE TimeInterval
#include <boost/test/unit_test.hpp>

#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(UnitScalarConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, 2.*rad/dm);
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    BOOST_TEST(
        interval.get_gradient_moment()
        == sycomore::Array<sycomore::Quantity>({20*rad/m,20*rad/m,20*rad/m}));
}

BOOST_AUTO_TEST_CASE(UnitVectorConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, {2*rad/dm,4*rad/m,8*rad/dam});
    BOOST_TEST(interval.get_duration() == 1.e-3*s);
    BOOST_TEST(
        interval.get_gradient_moment()
        == sycomore::Array<sycomore::Quantity>({20*rad/m,4*rad/m,0.8*rad/m}));
}
