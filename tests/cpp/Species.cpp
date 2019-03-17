#define BOOST_TEST_MODULE Species
#include <boost/test/unit_test.hpp>

#include "sycomore/Species.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(QuantityConstructorFull)
{
    using namespace sycomore::units;
    sycomore::Species const species(
        1.*ms, 10.*s, 3e-9*m*m/s, (1/1.5)/s, 0.9*rad/s, 1.1);
    BOOST_TEST(species.get_R1() == 1*kHz);
    BOOST_TEST(species.get_T1() == 1*ms);

    BOOST_TEST(species.get_R2() == 0.1*Hz);
    BOOST_TEST(species.get_T2() == 10*s);

    BOOST_TEST(species.get_D() == 3e-9*m*m/s);

    BOOST_TEST(species.get_R2_prime() == (1/1.5)/s);
    BOOST_TEST(species.get_T2_prime() == 1.5*s);

    BOOST_TEST(species.get_delta_omega() == 0.9*rad/s);

    BOOST_TEST(species.w == 1.1);
}

BOOST_AUTO_TEST_CASE(QuantityConstructorPartial)
{
    using namespace sycomore::units;
    sycomore::Species const species(1.*ms, 10.*s);
    BOOST_TEST(species.get_R1() == 1*kHz);
    BOOST_TEST(species.get_T1() == 1*ms);

    BOOST_TEST(species.get_R2() == 0.1*Hz);
    BOOST_TEST(species.get_T2() == 10*s);

    BOOST_TEST(species.get_D() == 0*m*m/s);

    BOOST_TEST(species.get_R2_prime() == 0*Hz);
//    BOOST_TEST(species.get_T2_prime() == INF*s);

    BOOST_TEST(species.get_delta_omega() == 0.*rad/s);

    BOOST_TEST(species.w == 1.);
}

