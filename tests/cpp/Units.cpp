#define BOOST_TEST_MODULE Units
#include <boost/test/unit_test.hpp>

#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(EmptyConstructor)
{
    sycomore::units::Time const duration;
    BOOST_TEST(duration.value == 0);
}

BOOST_AUTO_TEST_CASE(ValueConstructor)
{
    sycomore::units::Time const duration(1.);
    BOOST_TEST(duration.value == 1.);
}

BOOST_AUTO_TEST_CASE(Mult)
{
    using namespace sycomore::units;
    mult<Time, ElectricCurrent> x;
    Unit<0,0,1,1,0,0,0> y = x;
}

BOOST_AUTO_TEST_CASE(Div)
{
    using namespace sycomore::units;
    div<Length, Time> x;
    Unit<1,0,-1,0,0,0,0> y = x;
}

BOOST_AUTO_TEST_CASE(Pow)
{
    using namespace sycomore::units;
    // WARNING: unqualified "pow" is ambiguous, as math.h has been included,
    // possibly by boost.
    sycomore::units::pow<Length, 2> x;
    Unit<2,0,0,0,0,0,0> y = x;
}

BOOST_AUTO_TEST_CASE(Addition)
{
    using namespace sycomore::units;
    Length const one(1);
    Length const two(2);
    Length const three(one+two);
    BOOST_TEST(three.value == 3.);
}

BOOST_AUTO_TEST_CASE(Subtraction)
{
    using namespace sycomore::units;
    Length const one(1);
    Length const two(2);
    Length const minus_one(one-two);
    BOOST_TEST(minus_one.value == -1.);
}

BOOST_AUTO_TEST_CASE(ProductScalarLeft)
{
    using namespace sycomore::units;
    Length const two(2);
    Length const six(3*two);
    BOOST_TEST(six.value == 6.);
}

BOOST_AUTO_TEST_CASE(ProductScalarRight)
{
    using namespace sycomore::units;
    Length const two(2);
    Length const six(two*3);
    BOOST_TEST(six.value == 6.);
}

BOOST_AUTO_TEST_CASE(QuotientScalarLeft)
{
    using namespace sycomore::units;
    Time const six(6);
    Frequency const two = 12/six;
    BOOST_TEST(two.value == 2.);
}

BOOST_AUTO_TEST_CASE(ProductUnits)
{
    using namespace sycomore::units;
    Length const two(2);
    Length const three(3);
    Surface const six(two*three);
    BOOST_TEST(six.value == 6.);
}

BOOST_AUTO_TEST_CASE(QuotientUnits)
{
    using namespace sycomore::units;
    Length const six(6);
    Time const three(3);
    Velocity const two(six/three);
    BOOST_TEST(two.value == 2.);
}

BOOST_AUTO_TEST_CASE(BasicLiteral)
{
    using namespace sycomore::units;

    Length const length = 100_cm;
    BOOST_TEST(length.value == 1.);
}

BOOST_AUTO_TEST_CASE(DefinedBasicUnit)
{
    using namespace sycomore::units;

    Length const length = 100*cm;
    BOOST_TEST(length.value == 1.);
}

BOOST_AUTO_TEST_CASE(DerivedLiteral)
{
    using namespace sycomore::units;

    Force const force = 1_N;
    BOOST_TEST(force.value == 1.);
}

BOOST_AUTO_TEST_CASE(DefinedDerivedUnit)
{
    using namespace sycomore::units;

    Force const force = 1*kN;
    BOOST_TEST(force.value == 1000.);
}
