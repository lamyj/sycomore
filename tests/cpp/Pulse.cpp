#define BOOST_TEST_MODULE Pulse
#include <boost/test/unit_test.hpp>

#include "sycomore/Pulse.h"

BOOST_AUTO_TEST_CASE(Constructor)
{
    sycomore::Pulse const pulse(1, 2);
    BOOST_REQUIRE_EQUAL(pulse.angle, 1);
    BOOST_REQUIRE_EQUAL(pulse.phase, 2);
}

BOOST_AUTO_TEST_CASE(RotationMatrix)
{
    using namespace sycomore::units;
    sycomore::Pulse const pulse{41_deg, 27_deg};
    auto const m = pulse.rotation_matrix();

    // Values from CoMoTk@2b0ef02
    using sycomore::Complex;
    BOOST_TEST(std::norm(m[{0,0}] - Complex(8.773547901113861e-01, 0)) < 1e-9);
    BOOST_TEST(std::norm(m[{0,1}] - Complex(2.106079126622500e-01, -4.133413019334431e-01)) < 1e-9);
    BOOST_TEST(std::norm(m[{0,2}] - Complex(7.208904563684231e-02, 9.922205907857107e-02)) < 1e-9);

    BOOST_TEST(std::norm(m[{1,0}] - Complex(-2.106079126622500e-01, -4.133413019334431e-01)) < 1e-9);
    BOOST_TEST(std::norm(m[{1,1}] - Complex(7.547095802227720e-01, 0)) < 1e-9);
    BOOST_TEST(std::norm(m[{1,2}] - Complex(-2.106079126622500e-01, 4.133413019334431e-01)) < 1e-9);

    BOOST_TEST(std::norm(m[{2,0}] - Complex(7.208904563684231e-02, -9.922205907857107e-02)) < 1e-9);
    BOOST_TEST(std::norm(m[{2,1}] - Complex(2.106079126622500e-01, 4.133413019334431e-01)) < 1e-9);
    BOOST_TEST(std::norm(m[{2,2}] - Complex(8.773547901113861e-01, 0)) < 1e-9);
}
