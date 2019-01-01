#define BOOST_TEST_MODULE TimeInterval
#include <boost/test/unit_test.hpp>

#include "sycomore/TimeInterval.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(ScalarScalarConstructor)
{
    sycomore::TimeInterval const interval(1., 2.);
    BOOST_TEST(interval.duration == 1.);
    BOOST_TEST(interval.gradient_moment == std::vector<sycomore::Real>({2,2,2}));
}

BOOST_AUTO_TEST_CASE(UnitScalarConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, 2.);
    BOOST_TEST(interval.duration == 1.e-3);
    BOOST_TEST(interval.gradient_moment == std::vector<sycomore::Real>({2,2,2}));
}

BOOST_AUTO_TEST_CASE(ScalarVectorConstructor)
{
    sycomore::TimeInterval const interval(1., {2,3,4});
    BOOST_TEST(interval.duration == 1.);
    BOOST_TEST(interval.gradient_moment == std::vector<sycomore::Real>({2,3,4}));
}

BOOST_AUTO_TEST_CASE(UnitVectorConstructor)
{
    using namespace sycomore::units;
    sycomore::TimeInterval const interval(1._ms, {2,3,4});
    BOOST_TEST(interval.duration == 1.e-3);
    BOOST_TEST(interval.gradient_moment == std::vector<sycomore::Real>({2,3,4}));
}
