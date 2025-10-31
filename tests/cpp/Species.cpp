#define BOOST_TEST_MODULE Species
#include <boost/test/unit_test.hpp>

#include "sycomore/Species.h"
#include "sycomore/units.h"

#include "utils.h"

BOOST_AUTO_TEST_CASE(QuantityConstructorFull)
{
    using namespace sycomore::units;
    sycomore::Species const species(
        1.*ms, 10.*s, 3e-9*m*m/s, 0.9*rad/s);
    CHECK_QUANTITY(species.R1(), 1*kHz);
    CHECK_QUANTITY(species.T1(), 1*ms);

    CHECK_QUANTITY(species.R2(), 0.1*Hz);
    CHECK_QUANTITY(species.T2(), 10*s);

    CHECK_QUANTITY(species.D()[0], 3e-9*m*m/s);

    CHECK_QUANTITY(species.delta_omega(), 0.9*rad/s);
}

BOOST_AUTO_TEST_CASE(QuantityConstructorPartial)
{
    using namespace sycomore::units;
    sycomore::Species const species(1.*ms, 10.*s);
    CHECK_QUANTITY(species.R1(), 1*kHz);
    CHECK_QUANTITY(species.T1(), 1*ms);

    CHECK_QUANTITY(species.R2(), 0.1*Hz);
    CHECK_QUANTITY(species.T2(), 10*s);

    CHECK_QUANTITY(species.D()[0], 0*m*m/s);

    CHECK_QUANTITY(species.delta_omega(), 0.*rad/s);
}

BOOST_AUTO_TEST_CASE(DScalar)
{
    using namespace sycomore::units;
    sycomore::Species species(1.*ms, 10.*s);
    species.set_D(1*um*um/ms);
    sycomore::Matrix3x3Q const D{
        {1*um*um/ms, 0*um*um/ms, 0*um*um/ms},
        {0*um*um/ms, 1*um*um/ms, 0*um*um/ms},
        {0*um*um/ms, 0*um*um/ms, 1*um*um/ms}};
    CHECK_QUANTITY(species.D(), D);
}

BOOST_AUTO_TEST_CASE(DTensor)
{
    using namespace sycomore::units;
    sycomore::Species species(1.*ms, 10.*s);
    sycomore::Matrix3x3Q const D{
        {1*um*um/ms, 4*um*um/ms, 7*um*um/ms},
        {2*um*um/ms, 5*um*um/ms, 8*um*um/ms},
        {3*um*um/ms, 6*um*um/ms, 9*um*um/ms}};
    species.set_D(D);
    CHECK_QUANTITY(species.D(), D);
}
