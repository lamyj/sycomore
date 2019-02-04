#define BOOST_TEST_MODULE Species
#include <boost/test/unit_test.hpp>

#include "sycomore/Species.h"
#include "sycomore/units.h"

BOOST_AUTO_TEST_CASE(ScalarConstructorFull)
{
    using namespace sycomore::units;
    sycomore::Species const species(1., 10., 3e-9, 1/1.5, 0.9, 1.1);
    BOOST_TEST(species.R1 == 1);
    BOOST_TEST(species.R2 == 10);
    BOOST_TEST(species.D == 3e-9);
    BOOST_TEST(species.R2_prime == 1/1.5);
    BOOST_TEST(species.delta_omega == 0.9);
    BOOST_TEST(species.w == 1.1);
}

BOOST_AUTO_TEST_CASE(ScalarConstructorPartial)
{
    using namespace sycomore::units;
    sycomore::Species const species(1., 10.);
    BOOST_TEST(species.R1 == 1);
    BOOST_TEST(species.R2 == 10);
    BOOST_TEST(species.D == 0);
    BOOST_TEST(species.R2_prime == 0);
    BOOST_TEST(species.delta_omega == 0);
    BOOST_TEST(species.w == 1);
}

BOOST_AUTO_TEST_CASE(QuantityConstructorFull)
{
    using namespace sycomore::units;
    sycomore::Species const species(
        1000_ms, 1/100_ms, 3_um*um/ms, 1/1.5_s, 0.9_rad/s, 1.1);
    BOOST_TEST(species.R1 == 1);
    BOOST_TEST(species.R2 == 10);
    BOOST_TEST(species.D == 3e-9);
    BOOST_TEST(species.R2_prime == 1/1.5);
    BOOST_TEST(species.delta_omega == 0.9);
    BOOST_TEST(species.w == 1.1);
}

BOOST_AUTO_TEST_CASE(QuantityConstructorPartial)
{
    using namespace sycomore::units;
    sycomore::Species const species(1000_ms, 1/100_ms);
    BOOST_TEST(species.R1 == 1);
    BOOST_TEST(species.R2 == 10);
    BOOST_TEST(species.D == 0);
    BOOST_TEST(species.R2_prime == 0);
    BOOST_TEST(species.delta_omega == 0);
    BOOST_TEST(species.w == 1);
}
