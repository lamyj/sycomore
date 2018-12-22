#define BOOST_TEST_MODULE magnetization
#include <boost/test/unit_test.hpp>

#include "sycomore/magnetization.h"

BOOST_AUTO_TEST_CASE(ComplexMagnetizationFromMagnetization, *boost::unit_test::tolerance(1e-9))
{
    sycomore::Magnetization const m_r(1,2,3);
    auto && m_c = sycomore::as_complex_magnetization(m_r);

    BOOST_TEST(m_c(0).real() == std::sqrt(2.)/2.);
    BOOST_TEST(m_c(0).imag() == std::sqrt(2.));

    BOOST_TEST(m_c(1).real() == 3);
    BOOST_TEST(m_c(1).imag() == 0);

    BOOST_TEST(m_c(2).real() == std::sqrt(2.)/2.);
    BOOST_TEST(m_c(2).imag() == -std::sqrt(2.));
}

BOOST_AUTO_TEST_CASE(ComplexMagnetizationToMagnetization, *boost::unit_test::tolerance(1e-9))
{
    sycomore::ComplexMagnetization const m_c{
        {std::sqrt(2.)/2., std::sqrt(2.)},
        3,
        {std::sqrt(2.)/2., -std::sqrt(2.)}};

    auto && m_r = sycomore::as_real_magnetization(m_c);

    BOOST_TEST(m_r(0) == 1);
    BOOST_TEST(m_r(1) == 2);
    BOOST_TEST(m_r(2) == 3);
}
