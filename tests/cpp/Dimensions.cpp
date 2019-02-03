#define BOOST_TEST_MODULE Dimensions
#include <boost/test/unit_test.hpp>

#include "sycomore/Dimensions.h"

BOOST_AUTO_TEST_CASE(Comparison)
{
    sycomore::Dimensions const Length{1,0,0,0,0,0,0};
    sycomore::Dimensions const Time{0,0,1,0,0,0,0};
    BOOST_CHECK(sycomore::Length == sycomore::Length);
    BOOST_CHECK(sycomore::Length != sycomore::Time);
}

BOOST_AUTO_TEST_CASE(MultiplicationInPlace)
{
    sycomore::Dimensions d1{1,0,0,0,0,0,0};
    sycomore::Dimensions const d2{0,1,0,0,0,0,0};
    sycomore::Dimensions const r{1,1,0,0,0,0,0};
    d1 *= d2;
    BOOST_CHECK(d1 == r);
}

BOOST_AUTO_TEST_CASE(DivisionInPlace)
{
    sycomore::Dimensions d1{1,0,0,0,0,0,0};
    sycomore::Dimensions const d2{0,1,0,0,0,0,0};
    sycomore::Dimensions const r{1,-1,0,0,0,0,0};
    d1 /= d2;
    BOOST_CHECK(d1 == r);
}

BOOST_AUTO_TEST_CASE(Multiplication)
{
    sycomore::Dimensions const d1{1,0,0,0,0,0,0};
    sycomore::Dimensions const d2{0,1,0,0,0,0,0};
    sycomore::Dimensions const r{1,1,0,0,0,0,0};
    BOOST_CHECK(d1*d2 == r);
}

BOOST_AUTO_TEST_CASE(Division)
{
    sycomore::Dimensions const d1{1,0,0,0,0,0,0};
    sycomore::Dimensions const d2{0,1,0,0,0,0,0};
    sycomore::Dimensions const r{1,-1,0,0,0,0,0};
    BOOST_CHECK(d1/d2 == r);
}

BOOST_AUTO_TEST_CASE(Pow)
{
    sycomore::Dimensions const d{3,0,-2,0,0,0,0};
    sycomore::Dimensions const r{6,0,-4,0,0,0,0};
    BOOST_CHECK(std::pow(d, 2) == r);
}
