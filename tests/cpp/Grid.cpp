#define BOOST_TEST_MODULE Grid
#include <boost/test/unit_test.hpp>

#include "sycomore/Grid.h"
#include "sycomore/GridScanner.h"
#include "sycomore/magnetization.h"

using Grid = sycomore::Grid<sycomore::ComplexMagnetization>;

BOOST_AUTO_TEST_CASE(EmptyConstructor)
{
    Grid const grid;
    BOOST_TEST(grid.dimension() == 0);
    BOOST_TEST(grid.origin().empty());
    BOOST_TEST(grid.shape().empty());
    BOOST_TEST(grid.stride().empty());
}

BOOST_AUTO_TEST_CASE(UninitializedConstructor)
{
    Grid const grid({-2,-3,-5}, {5,7,11});
    BOOST_TEST(grid.dimension() == 3);
    BOOST_TEST(grid.origin() == sycomore::Index({-2,-3,-5}));
    BOOST_TEST(grid.shape() == sycomore::Shape({5,7,11}));
    BOOST_TEST(grid.stride() == sycomore::Stride({1,5,35,385}));
}

BOOST_AUTO_TEST_CASE(InitializedConstructor)
{
    sycomore::ComplexMagnetization const m{{1, 2}, 3, {4, 5}};
    Grid const grid({-5, -3}, {11, 7}, m);
    BOOST_TEST(grid.dimension() == 2);
    BOOST_TEST(grid.origin() == sycomore::Index({-5, -3}));
    BOOST_TEST(grid.shape() == sycomore::Shape({11,7}));
    BOOST_TEST(grid.stride() == sycomore::Stride({1,11,77}));
    for(int i=0; i<grid.stride()[grid.dimension()]; ++i)
    {
        BOOST_REQUIRE(grid.data()[i] == m);
    }
}

BOOST_AUTO_TEST_CASE(Accessor)
{
    sycomore::ComplexMagnetization const zero{0,0,0};
    sycomore::ComplexMagnetization const non_zero{1,2,3};
    Grid grid({-5, -3}, {11, 7}, zero);
    grid[{-2,-1}] = non_zero;

    auto && grid_const = const_cast<Grid const &>(grid);

    for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
    {
        if(index.first == sycomore::Index{-2,-1})
        {
            BOOST_REQUIRE(grid_const[index.first] == non_zero);
        }
        else
        {
            BOOST_REQUIRE(grid_const[index.first] == zero);
        }
    }
}

BOOST_AUTO_TEST_CASE(ScanOrder)
{
    Grid grid({-5, -3}, {11, 7}, {0, 0, 0});

    unsigned int i=0;
    for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
    {
        grid[index.first] = {
            sycomore::Real(i), sycomore::Real(i), sycomore::Real(i)};
        ++i;
    }

    auto && data = grid.data();
    for(unsigned int i=0; i<grid.stride()[grid.dimension()]; ++i)
    {
        BOOST_REQUIRE(*data == sycomore::ComplexMagnetization({
            sycomore::Real(i), sycomore::Real(i), sycomore::Real(i)}));
        ++data;
    }
}

void compare_grids(
    Grid const & old, Grid const & new_,
    bool is_initialized, Grid::value_type const & value={})
{
    sycomore::GridScanner const generator_old(old.origin(), old.shape());
    std::vector<sycomore::Index> indices_old;
    std::transform(
        generator_old.begin(), generator_old.end(),
        std::back_inserter(indices_old),
        [](sycomore::GridScanner::value_type const & x){ return x.first; });

    sycomore::GridScanner const generator_new(new_.origin(), new_.shape());
    std::vector<sycomore::Index> indices_new;
    std::transform(
        generator_new.begin(), generator_new.end(),
        std::back_inserter(indices_new),
        [](sycomore::GridScanner::value_type const & x){ return x.first; });

    std::vector<sycomore::Index> intersection;
    std::vector<sycomore::Index> difference;
    for(auto && index: indices_new)
    {
        auto && it = std::find(indices_old.begin(), indices_old.end(), index);
        if(it != indices_old.end())
        {
            intersection.push_back(index);
        }
        else
        {
            difference.push_back(index);
        }
    }

    for(auto && index: intersection)
    {
        BOOST_REQUIRE(old[index] == new_[index]);
    }
    if(is_initialized)
    {
        for(auto && index: difference)
        {
            BOOST_REQUIRE(new_[index] == value);
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeLargerUninitialized)
{
    Grid old({-2,-1}, {5,3});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({-3,-2}, {7,5});

    BOOST_TEST(new_.origin() == sycomore::Index({-3,-2}));
    BOOST_TEST(new_.shape() == sycomore::Shape({7, 5}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 7, 35}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeLargerInitialized)
{
    Grid old({-2,-1}, {5,3});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({-3,-2}, {7,5}, {1000, 0, 0});

    BOOST_TEST(new_.origin() == sycomore::Index({-3,-2}));
    BOOST_TEST(new_.shape() == sycomore::Shape({7, 5}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 7, 35}));
    compare_grids(old, new_, true, {1000,0,0});
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerUninitialized)
{
    Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({-2,-1}, {5,3});

    BOOST_TEST(new_.origin() == sycomore::Index({-2,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,3}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 5, 15}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerInitialized)
{
    Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({-2,-1}, {5,3}, {1000, 0, 0});

    BOOST_TEST(new_.origin() == sycomore::Index({-2,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,3}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 5, 15}));
    compare_grids(old, new_, true, {1000, 0, 0});
}

BOOST_AUTO_TEST_CASE(ReshapeMixedUninitialized)
{
    Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({-5,-1}, {11,3});

    BOOST_TEST(new_.origin() == sycomore::Index({-5,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({11,3}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 11, 33}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeMixedInitialized)
{
    Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({-5,-1}, {11,3}, {1000,0,0});

    BOOST_TEST(new_.origin() == sycomore::Index({-5,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({11,3}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 11, 33}));
    compare_grids(old, new_, true, {1000,0,0});
}

BOOST_AUTO_TEST_CASE(ReshapeDisjointUninitialized)
{
    Grid old({-3,-2}, {3, 2});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({1,9}, {5,7});

    BOOST_TEST(new_.origin() == sycomore::Index({1,9}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,7}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 5, 35}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeDisjointInitialized)
{
    Grid old({-3,-2}, {3, 2});

    int i=1;
    for(auto && index: sycomore::GridScanner(old.origin(), old.shape()))
    {
        old[index.first] = {i,0,0};
        ++i;
    }

    Grid new_(old);
    new_.reshape({1,9}, {5,7}, {1000,0,0});

    BOOST_TEST(new_.origin() == sycomore::Index({1,9}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,7}));
    BOOST_TEST(new_.stride() == sycomore::Stride({1, 5, 35}));
    compare_grids(old, new_, true, {1000,0,0});
}

BOOST_AUTO_TEST_CASE(ConstIterator)
{
    Grid grid({-3,-2}, {7,5});
    unsigned int i=1;
    for(auto && index: sycomore::GridScanner(grid.origin(), grid.shape()))
    {
        grid[index.first] = {i,0,0};
        ++i;
    }

    Grid const & grid_const = const_cast<Grid const &>(grid);

    i=1;
    for(auto && x: grid_const)
    {
        BOOST_REQUIRE(sycomore::ComplexMagnetization(i,0,0)==x);
        ++i;
    }
}

BOOST_AUTO_TEST_CASE(Iterator)
{
    Grid grid({-3,-2}, {7,5});
    unsigned int i=1;
    for(auto && x: grid)
    {
        x = {i,0,0};
        ++i;
    }

    Grid const & grid_const = const_cast<Grid const &>(grid);

    i=1;
    for(auto && x: grid_const)
    {
        BOOST_REQUIRE(sycomore::ComplexMagnetization(i,0,0)==x);
        ++i;
    }
}
