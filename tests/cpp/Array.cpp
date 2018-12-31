#define BOOST_TEST_MODULE Array
#include <boost/test/unit_test.hpp>

#include "sycomore/Array.h"

BOOST_AUTO_TEST_CASE(EmptyConstructor)
{
    sycomore::Array<float> const array;
    BOOST_TEST(array.dimension() == 0);
    BOOST_TEST(array.shape().empty());
    BOOST_TEST(array.stride().empty());
}

BOOST_AUTO_TEST_CASE(UninitializedConstructor)
{
    sycomore::Array<float> const array({2,3,5});
    BOOST_TEST(array.dimension() == 3);
    BOOST_TEST(array.shape() == sycomore::Shape({2,3,5}));
    BOOST_TEST(array.stride() == sycomore::Stride({1,2,6,30}));
}

BOOST_AUTO_TEST_CASE(InitializedConstructor)
{
    sycomore::Array<float> const array({5, 3}, 42);
    BOOST_TEST(array.dimension() == 2);
    BOOST_TEST(array.shape() == sycomore::Shape({5,3}));
    BOOST_TEST(array.stride() == sycomore::Stride({1,5,15}));
    for(int i=0; i<array.stride()[array.dimension()]; ++i)
    {
        BOOST_TEST(array.data()[i] == 42);
    }
}

BOOST_AUTO_TEST_CASE(Accessor)
{
    sycomore::Array<float> array({5, 3}, 0);
    array[{2,1}] = 42;

    sycomore::Array<float> const & array_const = const_cast<
            sycomore::Array<float> const &
        >(array);

    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        if(index == sycomore::Index{2,1})
        {
            BOOST_TEST(array_const[index] == 42);
        }
        else
        {
            BOOST_TEST(array_const[index] == 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(ScanOrder)
{
    sycomore::Array<unsigned int> array({3, 5, 7, 11}, 0);

    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        array[index] = i;
        ++i;
    }

    auto && data = array.data();
    for(unsigned int i=0; i<array.stride()[array.dimension()]; ++i)
    {
        BOOST_TEST(*data == i);
        ++data;
    }
}

template<typename TScalar>
void compare_arrays(
    sycomore::Array<TScalar> const & old_array,
    sycomore::Array<TScalar> const & new_array,
    bool is_initialized, TScalar const & value={})
{
    sycomore::IndexGenerator const generator_old(old_array.shape());
    std::vector<sycomore::Index> const indices_old(
        generator_old.begin(), generator_old.end());

    sycomore::IndexGenerator generator_new(new_array.shape());
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
        BOOST_TEST(old_array[index] == new_array[index]);
    }
    if(is_initialized)
    {
        for(auto && index: difference)
        {
            BOOST_TEST(new_array[index] == value);
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeLargerUninitialized)
{
    sycomore::Array<int> old_array({5, 3});
    unsigned int i=1;
    for(auto && index: sycomore::IndexGenerator(old_array.shape()))
    {
        old_array[index] = i;
        ++i;
    }

    auto new_array = old_array;
    new_array.reshape({7,5});

    BOOST_TEST(new_array.shape() == sycomore::Shape({7,5}));
    compare_arrays(old_array, new_array, false);
}

BOOST_AUTO_TEST_CASE(ReshapeLargerInitialized)
{
    sycomore::Array<int> old_array({5, 3});
    unsigned int i=1;
    for(auto && index: sycomore::IndexGenerator(old_array.shape()))
    {
        old_array[index] = i;
        ++i;
    }

    auto new_array = old_array;
    new_array.reshape({7,5}, 1000);

    BOOST_TEST(new_array.shape() == sycomore::Shape({7,5}));
    compare_arrays(old_array, new_array, true, 1000);
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerUninitialized)
{
    sycomore::Array<int> old_array({7,5});
    unsigned int i=1;
    for(auto && index: sycomore::IndexGenerator(old_array.shape()))
    {
        old_array[index] = i;
        ++i;
    }

    auto new_array = old_array;
    new_array.reshape({5,3});

    BOOST_TEST(new_array.shape() == sycomore::Shape({5, 3}));
    compare_arrays(old_array, new_array, false);
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerInitialized)
{
    sycomore::Array<int> old_array({7,5});
    unsigned int i=1;
    for(auto && index: sycomore::IndexGenerator(old_array.shape()))
    {
        old_array[index] = i;
        ++i;
    }

    auto new_array = old_array;
    new_array.reshape({5,3}, 1000);

    BOOST_TEST(new_array.shape() == sycomore::Shape({5, 3}));
    compare_arrays(old_array, new_array, true, 1000);
}

BOOST_AUTO_TEST_CASE(ReshapeMixedUninitialized)
{
    sycomore::Array<int> old_array({7,5});
    unsigned int i=1;
    for(auto && index: sycomore::IndexGenerator(old_array.shape()))
    {
        old_array[index] = i;
        ++i;
    }

    auto new_array = old_array;
    new_array.reshape({11,3});

    BOOST_TEST(new_array.shape() == sycomore::Shape({11, 3}));
    compare_arrays(old_array, new_array, false);
}

BOOST_AUTO_TEST_CASE(ReshapeMixedInitialized)
{
    sycomore::Array<int> old_array({7,5});
    unsigned int i=1;
    for(auto && index: sycomore::IndexGenerator(old_array.shape()))
    {
        old_array[index] = i;
        ++i;
    }

    auto new_array = old_array;
    new_array.reshape({11,3}, 1000);

    BOOST_TEST(new_array.shape() == sycomore::Shape({11, 3}));
    compare_arrays(old_array, new_array, true, 1000);
}
