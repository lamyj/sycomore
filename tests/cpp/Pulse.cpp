#define BOOST_TEST_MODULE Pulse
#include <boost/test/unit_test.hpp>

#include "sycomore/Pulse.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(Constructor)
{
    using namespace sycomore::units;
    sycomore::Pulse const pulse(1*rad, 2*rad);
    BOOST_REQUIRE_EQUAL(pulse.get_angle(), 1*rad);
    BOOST_REQUIRE_EQUAL(pulse.get_phase(), 2*rad);
}

BOOST_AUTO_TEST_CASE(RotationMatrix, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;
    sycomore::Pulse const pulse{41_deg, 27_deg};
    auto const m = pulse.rotation_matrix();

    // Values from CoMoTk@2b0ef02
    using sycomore::Complex;

    auto x = m[{0,0}];
    BOOST_TEST(x.real() == 8.773547901113861e-01);
    BOOST_TEST(x.imag() == 0);

    x = m[{0,1}];
    BOOST_TEST(x.real() == 2.106079126622500e-01);
    BOOST_TEST(x.imag() == -4.133413019334431e-01);

    x = m[{0,2}];
    BOOST_TEST(x.real() == 7.208904563684231e-02);
    BOOST_TEST(x.imag() == 9.922205907857107e-02);

    x = m[{1,0}];
    BOOST_TEST(x.real() = -2.106079126622500e-01);
    BOOST_TEST(x.imag() = -4.133413019334431e-01);

    x = m[{1,0}];
    BOOST_TEST(x.real() == -2.106079126622500e-01);
    BOOST_TEST(x.imag() == -4.133413019334431e-01);

    x = m[{1,1}];
    BOOST_TEST(x.real() == 7.547095802227720e-01);
    BOOST_TEST(x.imag() == 0.);

    x = m[{1,2}];
    BOOST_TEST(x.real() == -2.106079126622500e-01);
    BOOST_TEST(x.imag() == 4.133413019334431e-01);

    x = m[{2,0}];
    BOOST_TEST(x.real() == 7.208904563684231e-02);
    BOOST_TEST(x.imag() == -9.922205907857107e-02);

    x = m[{2,1}];
    BOOST_TEST(x.real() == 2.106079126622500e-01);
    BOOST_TEST(x.imag() == 4.133413019334431e-01);

    x = m[{2,2}];
    BOOST_TEST(x.real() == 8.773547901113861e-01);
    BOOST_TEST(x.imag() == 0);
}
