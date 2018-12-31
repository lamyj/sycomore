#define BOOST_TEST_MODULE Grid
#include <boost/test/unit_test.hpp>

#include "sycomore/Grid.h"

BOOST_AUTO_TEST_CASE(EmptyConstructor)
{
    sycomore::Grid const grid;
    BOOST_TEST(grid.dimension() == 0);
    BOOST_TEST(grid.origin().empty());
    BOOST_TEST(grid.shape().empty());
    BOOST_TEST(grid.stride().empty());
}

BOOST_AUTO_TEST_CASE(UninitializedConstructor)
{
    sycomore::Grid const grid({-2,-3,-5}, {5,7,11});
    BOOST_TEST(grid.dimension() == 3);
    BOOST_TEST(grid.origin() == sycomore::Index({-2,-3,-5}));
    BOOST_TEST(grid.shape() == sycomore::Shape({5,7,11}));
    BOOST_TEST(grid.stride() == sycomore::Stride({1,5,35,385}));
}

BOOST_AUTO_TEST_CASE(InitializedConstructor)
{
    sycomore::ComplexMagnetization const m{{1, 2}, 3, {4, 5}};
    sycomore::Grid const grid({-5, -3}, {11, 7}, m);
    BOOST_TEST(grid.dimension() == 2);
    BOOST_TEST(grid.origin() == sycomore::Index({-5, -3}));
    BOOST_TEST(grid.shape() == sycomore::Shape({11,7}));
    BOOST_TEST(grid.stride() == sycomore::Stride({1,11,77}));
    for(int i=0; i<grid.stride()[grid.dimension()]; ++i)
    {
        BOOST_TEST(grid.data()[i] == m);
    }
}

BOOST_AUTO_TEST_CASE(Accessor)
{
    sycomore::ComplexMagnetization const zero{0,0,0};
    sycomore::ComplexMagnetization const non_zero{1,2,3};
    sycomore::Grid grid({-5, -3}, {11, 7}, zero);
    grid[{-2,-1}] = non_zero;

    auto && grid_const = const_cast<sycomore::Grid const &>(grid);

    for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
    {
        if(index == sycomore::Index{-2,-1})
        {
            BOOST_TEST(grid_const[index] == non_zero);
        }
        else
        {
            BOOST_TEST(grid_const[index] == zero);
        }
    }
}

BOOST_AUTO_TEST_CASE(ScanOrder)
{
    sycomore::Grid grid({-5, -3}, {11, 7}, {0, 0, 0});

    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(grid.origin(), grid.shape()))
    {
        grid[index] = {i, i, i};
        ++i;
    }

    auto && data = grid.data();
    for(unsigned int i=0; i<grid.stride()[grid.dimension()]; ++i)
    {
        BOOST_TEST(*data == sycomore::ComplexMagnetization({i,i,i}));
        ++data;
    }
}

void compare_grids(
    sycomore::Grid const & old, sycomore::Grid const & new_,
    bool is_initialized, sycomore::Grid::value_type const & value={})
{
    sycomore::IndexGenerator const generator_old(old.origin(), old.shape());
    std::vector<sycomore::Index> const indices_old(
        generator_old.begin(), generator_old.end());

    sycomore::IndexGenerator const generator_new(new_.origin(), new_.shape());
    std::vector<sycomore::Index> const indices_new(
        generator_new.begin(), generator_new.end());

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
        BOOST_TEST(old[index] == new_[index]);
    }
    if(is_initialized)
    {
        for(auto && index: difference)
        {
            BOOST_TEST(new_[index] == value);
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeLargerUninitialized)
{
    sycomore::Grid old({-2,-1}, {5,3});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({-3,-2}, {7,5});

    BOOST_TEST(new_.origin() == sycomore::Index({-3,-2}));
    BOOST_TEST(new_.shape() == sycomore::Shape({7, 5}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeLargerInitialized)
{
    sycomore::Grid old({-2,-1}, {5,3});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({-3,-2}, {7,5}, {1000, 0, 0});

    BOOST_TEST(new_.origin() == sycomore::Index({-3,-2}));
    BOOST_TEST(new_.shape() == sycomore::Shape({7, 5}));
    compare_grids(old, new_, true, {1000,0,0});
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerUninitialized)
{
    sycomore::Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({-2,-1}, {5,3});

    BOOST_TEST(new_.origin() == sycomore::Index({-2,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,3}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerInitialized)
{
    sycomore::Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({-2,-1}, {5,3}, {1000, 0, 0});

    BOOST_TEST(new_.origin() == sycomore::Index({-2,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,3}));
    compare_grids(old, new_, true, {1000, 0, 0});
}

BOOST_AUTO_TEST_CASE(ReshapeMixedUninitialized)
{
    sycomore::Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({-5,-1}, {11,3});

    BOOST_TEST(new_.origin() == sycomore::Index({-5,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({11,3}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeMixedInitialized)
{
    sycomore::Grid old({-3,-2}, {7, 5});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({-5,-1}, {11,3}, {1000,0,0});

    BOOST_TEST(new_.origin() == sycomore::Index({-5,-1}));
    BOOST_TEST(new_.shape() == sycomore::Shape({11,3}));
    compare_grids(old, new_, true, {1000,0,0});
}

BOOST_AUTO_TEST_CASE(ReshapeDisjointUninitialized)
{
    sycomore::Grid old({-3,-2}, {3, 2});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({1,9}, {5,7});

    BOOST_TEST(new_.origin() == sycomore::Index({1,9}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,7}));
    compare_grids(old, new_, false);
}

BOOST_AUTO_TEST_CASE(ReshapeDisjointInitialized)
{
    sycomore::Grid old({-3,-2}, {3, 2});

    int i=1;
    for(auto && index: sycomore::IndexGenerator(old.origin(), old.shape()))
    {
        old[index] = {i,0,0};
        ++i;
    }

    sycomore::Grid new_(old);
    new_.reshape({1,9}, {5,7}, {1000,0,0});

    BOOST_TEST(new_.origin() == sycomore::Index({1,9}));
    BOOST_TEST(new_.shape() == sycomore::Shape({5,7}));
    compare_grids(old, new_, true, {1000,0,0});
}
