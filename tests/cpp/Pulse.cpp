#define BOOST_TEST_MODULE Pulse
#include <boost/test/unit_test.hpp>

#include "sycomore/Pulse.h"

BOOST_AUTO_TEST_CASE(Constructor)
{
    sycomore::Pulse const pulse{1, 2};
    BOOST_REQUIRE_EQUAL(pulse.angle, 1);
    BOOST_REQUIRE_EQUAL(pulse.phase, 2);
}

BOOST_AUTO_TEST_CASE(RotationMatrix, *boost::unit_test::tolerance(1e-9))
{
    sycomore::Pulse const pulse{sycomore::deg2rad(41), sycomore::deg2rad(27)};
    auto const m = pulse.rotation_matrix();

    // Values from CoMoTk@2b0ef02
    BOOST_TEST(m(0,0).real() ==  8.773547901113861e-01); BOOST_TEST(m(0,0).imag() ==  0);
    BOOST_TEST(m(0,1).real() ==  2.106079126622500e-01); BOOST_TEST(m(0,1).imag() == -4.133413019334431e-01);
    BOOST_TEST(m(0,2).real() ==  7.208904563684231e-02); BOOST_TEST(m(0,2).imag() ==  9.922205907857107e-02);

    BOOST_TEST(m(1,0).real() == -2.106079126622500e-01); BOOST_TEST(m(1,0).imag() == -4.133413019334431e-01);
    BOOST_TEST(m(1,1).real() ==  7.547095802227720e-01); BOOST_TEST(m(1,1).imag() ==  0);
    BOOST_TEST(m(1,2).real() == -2.106079126622500e-01); BOOST_TEST(m(1,2).imag() ==  4.133413019334431e-01);

    BOOST_TEST(m(2,0).real() ==  7.208904563684231e-02); BOOST_TEST(m(2,0).imag() == -9.922205907857107e-02);
    BOOST_TEST(m(2,1).real() ==  2.106079126622500e-01); BOOST_TEST(m(2,1).imag() ==  4.133413019334431e-01);
    BOOST_TEST(m(2,2).real() ==  8.773547901113861e-01); BOOST_TEST(m(2,2).imag() ==  0);
}
