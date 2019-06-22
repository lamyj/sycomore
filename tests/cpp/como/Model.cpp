#define BOOST_TEST_MODULE como_Model
#include <boost/test/unit_test.hpp>

#include "sycomore/como/Model.h"
#include "sycomore/GridScanner.h"
#include "sycomore/Species.h"
#include "sycomore/TimeInterval.h"

void test_magnetization(
    sycomore::ComplexMagnetization const & m1,
    sycomore::ComplexMagnetization const & m2, sycomore::Real tolerance=1e-9)
{
    BOOST_TEST(m1.p.real() == m2.p.real(), boost::test_tools::tolerance(tolerance));
    BOOST_TEST(m1.p.imag() == m2.p.imag(), boost::test_tools::tolerance(tolerance));

    BOOST_TEST(m1.z == m2.z, boost::test_tools::tolerance(tolerance));

    BOOST_TEST(m1.m.real() == m2.m.real(), boost::test_tools::tolerance(tolerance));
    BOOST_TEST(m1.m.imag() == m2.m.imag(), boost::test_tools::tolerance(tolerance));
}

BOOST_AUTO_TEST_CASE(Constructor)
{
    using namespace sycomore::units;

    sycomore::Species const species(1_s, 0.1_s);
    sycomore::TimeInterval const echo(10_ms);

    sycomore::como::Model const model(species, {0,0,1}, {{"echo", echo}});

    BOOST_TEST((
        model.dimensions() == std::map<std::string, size_t>{{"echo", 0}}));
    BOOST_TEST((
        model.time_intervals() == std::vector<sycomore::TimeInterval>{{echo}}));

    BOOST_TEST((model.magnetization().origin() <= sycomore::Index{0}));
    BOOST_TEST((model.magnetization().shape() >= sycomore::Shape{1}));
}

BOOST_AUTO_TEST_CASE(Pulse)
{
    using namespace sycomore::units;

    sycomore::Species const species(1_s, 0.1_s);
    sycomore::como::Model model(species, {0,0,1}, {{"dummy", {0*s}}});

    sycomore::Pulse const pulse{41_deg, 27_deg};
    model.apply_pulse(pulse);

    auto && grid = model.magnetization();

    BOOST_TEST((grid.origin() <= sycomore::Index{0}));
    BOOST_TEST((grid.shape() >= sycomore::Shape{1}));

    for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
    {
        auto & m = grid[index.first];
        if(index.first == sycomore::Index{0})
        {
            // Values from CoMoTk@2b0ef02
            test_magnetization(
                m,
                {
                    {0.210607912662250, -0.413341301933443},
                    0.754709580222772,
                    {0.210607912662250, +0.413341301933443}
                });
        }
        else
        {
            test_magnetization(m, sycomore::ComplexMagnetization::zero);
        }
    }
}

BOOST_AUTO_TEST_CASE(TimeInterval, *boost::unit_test::tolerance(1e-9))
{
    using namespace sycomore::units;

    sycomore::Species const species{std::log(2)*Hz, std::log(2)*Hz};
    // Dummy time intervals to build a 2D model
    sycomore::como::Model model(species, {0, 0, 1}, {{"foo", {1*s}}, {"bar", {1*s}}});

    // Resulting complex magnetization: 1/2, sqrt(2)/2, 1/2
    model.apply_pulse({45_deg, 90_deg});

    // Transverse magnetization is divided by 2, and becomes 1/4; longitudinal
    // magnetization becomes 1/2*sqrt(2)/2 + (1-1/2) == 1/2*(1+sqrt(2)/2)
    // This value includes relaxation and repolarization effects with
    // m_eq == (0,0,1)
    model.apply_time_interval("foo");

    {
        auto && grid = model.magnetization();

        BOOST_TEST((model.magnetization().origin() <= sycomore::Index{-1, 0}));
        BOOST_TEST((model.magnetization().shape() >= sycomore::Shape{3, 0}));

        for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
        {
            auto && m = grid[index.first];
            if(index.first == sycomore::Index{-1,0})
            {
                test_magnetization(m, {0., 0., 0.25});
            }
            else if(index.first == sycomore::Index{0,0})
            {
                test_magnetization(m, {0., 1./2.*(1+std::sqrt(2.)/2.), 0.});
            }
            else if(index.first == sycomore::Index{1,0})
            {
                test_magnetization(m, {0.25, 0., 0.});
            }
            else
            {
                test_magnetization(m, sycomore::ComplexMagnetization::zero);
            }
        }
    }

    // Same modifications as before: transverse magnetization becomes 1/8,
    // longitudinal magnetization becomes 1/2 + 1/4*(1+sqrt(2)/2)
    // Configurations (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), and (0, 1) are
    // unoccupied
    model.apply_time_interval("bar");

    {
        auto && grid = model.magnetization();

        BOOST_TEST((model.magnetization().origin() <= sycomore::Index{-1, -1}));
        BOOST_TEST((model.magnetization().shape() >= sycomore::Shape{3, 3}));

        for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
        {
            auto && m = grid[index.first];
            if(index.first == sycomore::Index{-1,-1})
            {
                test_magnetization(m, {0., 0., 0.125});
            }
            else if(index.first == sycomore::Index{0,0})
            {
                test_magnetization(m, {0., 0.5 + 0.25*(1+std::sqrt(2.)/2.), 0.});
            }
            else if(index.first == sycomore::Index{1,1})
            {
                test_magnetization(m, {0.125, 0., 0.});
            }
            else
            {
                test_magnetization(m, sycomore::ComplexMagnetization::zero);
            }
        }
    }

    auto isochromat = model.isochromat();
    BOOST_TEST(isochromat[0] == 0.125*std::sqrt(2.));
    BOOST_TEST(isochromat[1] == 0.);
    BOOST_TEST(isochromat[2] == 0.5 + 0.25*(1+std::sqrt(2.)/2.));

    isochromat = model.isochromat({{0,0}, {-1,-1}});
    BOOST_TEST(isochromat[0] == 0.125*std::sqrt(2.)/2.);
    BOOST_TEST(isochromat[1] == 0.);
    BOOST_TEST(isochromat[2] == 0.5 + 0.25*(1+std::sqrt(2.)/2.));
}

BOOST_AUTO_TEST_CASE(Diffusion)
{
    using namespace sycomore::units;

    sycomore::Species const species{0_Hz, 0_Hz, 1_um*um/ms};
    sycomore::TimeInterval const echo{500_ms, 0.1*rad/um};

    sycomore::como::Model model(species, {0,0,1}, {{"echo", echo}});

    model.apply_pulse({40_deg, 0_deg});
    model.apply_time_interval("echo");

    // model = CoMoTk;
    // model.R1 = 0;
    // model.R2 = 0;
    // model.D = 1;
    // model.init_configuration([0;0;1]);
    // model.RF(deg2rad(40), 0);
    // model.time(1, 'tau', 500, 'p', 0.1);
    {
        auto && grid = model.magnetization();

        BOOST_TEST((model.magnetization().origin() <= sycomore::Index{-1}));
        BOOST_TEST((model.magnetization().shape() >= sycomore::Shape{3}));

        for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
        {
            auto && m = grid[index.first];
            if(index.first == sycomore::Index{-1})
            {
                test_magnetization(m, {0., 0., {0, 0.003062528150606}});
            }
            else if(index.first == sycomore::Index{0})
            {
                test_magnetization(m, {0., 0.766044443118978, 0.});
            }
            else if(index.first == sycomore::Index{1})
            {
                test_magnetization(m, {{0., -0.003062528150606}, 0., 0.});
            }
            else
            {
                test_magnetization(m, sycomore::ComplexMagnetization::zero);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(CleanUp)
{
    using namespace sycomore::units;

    sycomore::Species const species{std::log(2)*Hz, std::log(2)*Hz};
    // Dummy time intervals to build a 2D model
    sycomore::como::Model model(
        species, {0, 0, 1}, {{"short", {1*s}}, {"long", {2*s}}});
    model.set_epsilon(1.5*1./(1<<4));

    // Resulting complex magnetization: 1/2, sqrt(2)/2, 1/2
    model.apply_pulse({45_deg, 90_deg});
    model.apply_time_interval("long");

    // All occupied configurations are above the threshold
    {
        auto && grid = model.magnetization();

        BOOST_TEST((model.magnetization().origin() <= sycomore::Index{0, -1}));
        BOOST_TEST((model.magnetization().shape() >= sycomore::Shape{0, 1}));

        for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
        {
            auto && m = grid[index.first];
            if(index.first == sycomore::Index{0,-1})
            {
                test_magnetization(m, {0., 0., 0.125});
            }
            else if(index.first == sycomore::Index{0,0})
            {
                test_magnetization(m, {0., 0.25*std::sqrt(2.)/2. + 1-0.25, 0.});
            }
            else if(index.first == sycomore::Index{0,1})
            {
                test_magnetization(m, {0.125, 0., 0.});
            }
            else
            {
                test_magnetization(m, sycomore::ComplexMagnetization::zero);
            }
        }
    }

    model.apply_time_interval("short");

    // Magnetizations on non-zero configurations fall below the threshold
    {
        auto && grid = model.magnetization();

        BOOST_TEST((model.magnetization().origin() <= sycomore::Index{0, -3}));
        BOOST_TEST((model.magnetization().shape() >= sycomore::Shape{0, 3}));

        for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
        {
            auto && m = grid[index.first];
            if(index.first == sycomore::Index{0,0})
            {
                test_magnetization(
                    m, {0., 0.5*(0.25*std::sqrt(2.)/2. + 1-0.25)+(1-0.5), 0.});
            }
            else
            {
                test_magnetization(m, sycomore::ComplexMagnetization::zero);
            }
        }
    }
}
