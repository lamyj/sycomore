#define BOOST_TEST_MODULE Model
#include <boost/test/unit_test.hpp>

#include "sycomore/Model.h"
#include "sycomore/Species.h"
#include "sycomore/TimeInterval.h"

void test_magnetization(
    sycomore::ComplexMagnetization const & m1,
    sycomore::ComplexMagnetization const & m2, sycomore::Real tolerance=1e-9)
{
    BOOST_TEST(std::abs(m1.plus.real()-m2.plus.real())<=tolerance);
    BOOST_TEST(std::abs(m1.plus.imag()-m2.plus.imag())<=tolerance);

    BOOST_TEST(std::abs(m1.zero-m2.zero)<=tolerance);

    BOOST_TEST(std::abs(m1.minus.real()-m2.minus.real())<=tolerance);
    BOOST_TEST(std::abs(m1.minus.imag()-m2.minus.imag())<=tolerance);
}

BOOST_AUTO_TEST_CASE(Constructor)
{
    sycomore::Species const species(1*sycomore::units::s, 0.1*sycomore::units::s);
    sycomore::TimeInterval const echo(10*sycomore::units::ms);

    sycomore::Model const model(species, {0,0,1}, {{"echo", echo}});

    BOOST_TEST((model.dimensions() == std::map<std::string, size_t>{{"echo", 0}}));
    BOOST_TEST((
        model.time_intervals()
        == std::map<std::string, sycomore::TimeInterval>{{"echo", echo}}));

    BOOST_TEST((model.grid().origin() <= sycomore::Index{0}));
    BOOST_TEST((model.grid().shape() >= sycomore::Shape{1}));
}

BOOST_AUTO_TEST_CASE(Pulse)
{
    sycomore::Species const species(1*sycomore::units::s, 0.1*sycomore::units::s);
    sycomore::Model model(species, {0,0,1}, {{"dummy", {0}}});

    sycomore::Pulse const pulse{sycomore::deg2rad(30), sycomore::deg2rad(180)};
    model.apply_pulse(pulse);

    auto && grid = model.grid();

    BOOST_TEST((grid.origin() <= sycomore::Index{0}));
    BOOST_TEST((grid.shape() >= sycomore::Shape{1}));

    for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
    {
        auto & m = grid[index];
        if(index == sycomore::Index{0})
        {
            // Values from CoMoTk@2b0ef02
            test_magnetization(
                m,
                {{0., std::sqrt(2.)/4.}, std::sqrt(3.)/2., {0., -std::sqrt(2.)/4.}});
        }
        else
        {
            test_magnetization(m, {0., 0., 0.});
        }
    }
}

//for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
//{
//    std::cout << "("; for(auto && x: index) { std::cout << " " << x; }; std::cout << " ): ";
//    std::cout << grid[index].plus << " " << grid[index].zero << " " << grid[index].minus << std::endl;
//}

BOOST_AUTO_TEST_CASE(TimeInterval, *boost::unit_test::tolerance(1e-9))
{
    sycomore::Species const species{
        std::log(2)*sycomore::units::Hz, std::log(2)*sycomore::units::Hz};
    // Dummy time intervals to build a 2D model
    sycomore::Model model(species, {0, 0, 1}, {{"foo", {1}}, {"bar", {1}}});

    // Resulting complex magnetization: 1/2, sqrt(2)/2, 1/2
    model.apply_pulse({sycomore::deg2rad(45), sycomore::deg2rad(90)});

    // Transverse magnetization is divided by 2, and becomes 1/4; longitudinal
    // magnetization becomes 1/2*sqrt(2)/2 + (1-1/2) == 1/2*(1+sqrt(2)/2)
    // This value includes relaxation and repolarization effects with
    // m_eq == (0,0,1)
    model.apply_time_interval("foo");

    {
        auto && grid = model.grid();

        BOOST_TEST((model.grid().origin() <= sycomore::Index{-1, 0}));
        BOOST_TEST((model.grid().shape() >= sycomore::Shape{3, 0}));

        for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
        {
            auto && m = grid[index];
            if(index == sycomore::Index{-1,0})
            {
                test_magnetization(m, {0., 0., 0.25});
            }
            else if(index == sycomore::Index{0,0})
            {
                test_magnetization(m, {0., 1./2.*(1+std::sqrt(2.)/2.), 0.});
            }
            else if(index == sycomore::Index{1,0})
            {
                test_magnetization(m, {0.25, 0., 0.});
            }
            else
            {
                test_magnetization(m, {0., 0., 0.});
            }
        }
    }

    // Same modifications as before: transverse magnetization becomes 1/8,
    // longitudinal magnetization becomes 1/2 + 1/4*(1+sqrt(2)/2)
    // Configurations (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), and (0, 1) are
    // unoccupied
    model.apply_time_interval("bar");

    {
        auto && grid = model.grid();

        BOOST_TEST((model.grid().origin() <= sycomore::Index{-3, 0}));
        BOOST_TEST((model.grid().shape() >= sycomore::Shape{7, 0}));

        for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
        {
            auto && m = grid[index];
            if(index == sycomore::Index{-1,-1})
            {
                test_magnetization(m, {0., 0., 0.125});
            }
            else if(index == sycomore::Index{0,0})
            {
                test_magnetization(m, {0., 0.5 + 0.25*(1+std::sqrt(2.)/2.), 0.});
            }
            else if(index == sycomore::Index{1,1})
            {
                test_magnetization(m, {0.125, 0., 0.});
            }
            else
            {
                test_magnetization(m, {0., 0., 0.});
            }
        }
    }

    auto isochromat = model.isochromat();
    BOOST_TEST(isochromat.x == 0.125*std::sqrt(2.));
    BOOST_TEST(isochromat.y == 0.);
    BOOST_TEST(isochromat.z == 0.5 + 0.25*(1+std::sqrt(2.)/2.));

    isochromat = model.isochromat({{0,0}, {-2,-2}});
    BOOST_TEST(isochromat.x == 0.);
    BOOST_TEST(isochromat.y == 0.);
    BOOST_TEST(isochromat.z == 0.5 + 0.25*(1+std::sqrt(2.)/2.));
}

BOOST_AUTO_TEST_CASE(CleanUp)
{
    sycomore::Species const species{
        std::log(2)*sycomore::units::Hz, std::log(2)*sycomore::units::Hz};
    // Dummy time intervals to build a 2D model
    sycomore::Model model(species, {0, 0, 1}, {{"short", {1}}, {"long", {2}}});
    model.epsilon(1.5*1./(1<<4));

    model.apply_pulse({sycomore::deg2rad(45), sycomore::deg2rad(90)});
    model.apply_time_interval("long");

    // All occupied configurations are above the threshold
    {
        auto && grid = model.grid();

        BOOST_TEST((model.grid().origin() <= sycomore::Index{0, -1}));
        BOOST_TEST((model.grid().shape() >= sycomore::Shape{0, 1}));

        for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
        {
            auto && m = grid[index];
            if(index == sycomore::Index{0,-1})
            {
                test_magnetization(m, {0., 0., 0.125});
            }
            else if(index == sycomore::Index{0,0})
            {
                test_magnetization(m, {0., 0.25*std::sqrt(2.)/2. + 1-0.25, 0.});
            }
            else if(index == sycomore::Index{0,1})
            {
                test_magnetization(m, {0.125, 0., 0.});
            }
            else
            {
                test_magnetization(m, {0., 0., 0.});
            }
        }
    }

    model.apply_time_interval("short");

    // Magnetizations on non-zero configurations fall below the threshold
    {
        auto && grid = model.grid();

        BOOST_TEST((model.grid().origin() <= sycomore::Index{0, -3}));
        BOOST_TEST((model.grid().shape() >= sycomore::Shape{0, 3}));

        for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
        {
            auto && m = grid[index];
            if(index == sycomore::Index{0,0})
            {
                test_magnetization(
                    m, {0., 0.5*(0.25*std::sqrt(2.)/2. + 1-0.25)+(1-0.5), 0.});
            }
            else
            {
                test_magnetization(m, {0., 0., 0.});
            }
        }
    }
}
