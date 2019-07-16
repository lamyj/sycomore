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

    BOOST_TEST(species.get_D()[0] == 3e-9*m*m/s);

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

    BOOST_TEST(species.get_D()[0] == 0*m*m/s);

    BOOST_TEST(species.get_R2_prime() == 0*Hz);
//    BOOST_TEST(species.get_T2_prime() == INF*s);

    BOOST_TEST(species.get_delta_omega() == 0.*rad/s);

    BOOST_TEST(species.w == 1.);
}

BOOST_AUTO_TEST_CASE(DScalar)
{
    using namespace sycomore::units;
    sycomore::Species species(1.*ms, 10.*s);
    species.set_D(1*um*um/ms);
    sycomore::Array<sycomore::Quantity> const D{
        1*um*um/ms, 0*um*um/ms, 0*um*um/ms,
        0*um*um/ms, 1*um*um/ms, 0*um*um/ms,
        0*um*um/ms, 0*um*um/ms, 1*um*um/ms};
    BOOST_TEST(species.get_D() == D);
}

BOOST_AUTO_TEST_CASE(DTensor)
{
    using namespace sycomore::units;
    sycomore::Species species(1.*ms, 10.*s);
    sycomore::Array<sycomore::Quantity> const D{
        1*um*um/ms, 4*um*um/ms, 7*um*um/ms,
        2*um*um/ms, 5*um*um/ms, 8*um*um/ms,
        3*um*um/ms, 6*um*um/ms, 9*um*um/ms};
    species.set_D(D);
    BOOST_TEST(species.get_D() == D);
}
