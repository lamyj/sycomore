#define BOOST_TEST_MODULE Units
#include <boost/test/unit_test.hpp>

#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(BasicLiteral)
{
    using namespace sycomore::units;

    sycomore::Quantity const length = 100_cm;
    BOOST_TEST(length.magnitude == 1.);
    BOOST_TEST(length.dimensions == sycomore::Length);
}

BOOST_AUTO_TEST_CASE(DefinedBasicUnit)
{
    using namespace sycomore::units;

    sycomore::Quantity const length = 100*cm;
    BOOST_TEST(length.magnitude == 1.);
    BOOST_TEST(length.dimensions == sycomore::Length);
}

BOOST_AUTO_TEST_CASE(DerivedLiteral)
{
    using namespace sycomore::units;

    sycomore::Quantity const force = 1_N;
    BOOST_TEST(force.magnitude == 1.);
    BOOST_TEST(force.dimensions == sycomore::Force);
}

BOOST_AUTO_TEST_CASE(DefinedDerivedUnit)
{
    using namespace sycomore::units;

    sycomore::Quantity const force = 1*kN;
    BOOST_TEST(force.magnitude == 1000.);
    BOOST_TEST(force.dimensions == sycomore::Force);
}
