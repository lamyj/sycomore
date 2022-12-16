#define BOOST_TEST_MODULE Pulse
#include <boost/test/unit_test.hpp>

#include "sycomore/Pulse.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(Constructor)
{
    using namespace sycomore::units;
    sycomore::Pulse const pulse_1(1*rad, 2*rad);
    BOOST_REQUIRE_EQUAL(pulse_1.get_angle(), 1*rad);
    BOOST_REQUIRE_EQUAL(pulse_1.get_phase(), 2*rad);
    
    sycomore::Pulse const pulse_2(2*rad);
    BOOST_REQUIRE_EQUAL(pulse_2.get_angle(), 2*rad);
    BOOST_REQUIRE_EQUAL(pulse_2.get_phase(), 0*rad);
}
