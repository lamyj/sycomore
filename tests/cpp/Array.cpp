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
    for(int i=0; i<array.stride()[array.stride().size()-1]; ++i)
    {
        BOOST_TEST(array.data()[i] == 42);
    }
}

BOOST_AUTO_TEST_CASE(IndexGenerator)
{
    sycomore::IndexGenerator generator({3,5});
    auto && iterator = generator.begin();

    for(unsigned int y=0; y<5; ++y)
    {
        for(unsigned int x=0; x<3; ++x)
        {
            BOOST_CHECK(iterator != generator.end());
            auto && index = *iterator;
            BOOST_CHECK(index == sycomore::Index({x,y}));
            ++iterator;
        }
    }

    BOOST_CHECK(iterator == generator.end());
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

BOOST_AUTO_TEST_CASE(ReshapeLargerUninitialized)
{
    sycomore::Shape const old_shape = {5, 3};
    sycomore::Shape const new_shape = {7, 5};

    sycomore::Array<float> array(old_shape);
    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        array[index] = i;
        ++i;
    }

    array.reshape(new_shape);

    BOOST_TEST(array.shape() == new_shape);

    sycomore::Index index(array.dimension(), 0);
    i=0;
    for(index[1]=0; index[1]<old_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<old_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == i);
            ++i;
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeLargerInitialized)
{
    sycomore::Shape const old_shape = {5, 3};
    sycomore::Shape const new_shape = {7, 5};

    sycomore::Array<float> array(old_shape);
    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        array[index] = i;
        ++i;
    }

    array.reshape(new_shape, 1000);

    BOOST_TEST(array.shape() == new_shape);

    sycomore::Index index(array.dimension(), 0);
    i=0;
    for(index[1]=0; index[1]<old_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<old_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == i);
            ++i;
        }
        for(index[0]=old_shape[0]; index[0]<new_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == 1000);
        }
    }
    for(index[1]=old_shape[1]; index[1]<new_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<new_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == 1000);
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerUninitialized)
{
    sycomore::Shape const old_shape{7, 5};
    sycomore::Shape const new_shape{5, 3};

    sycomore::Array<float> array(old_shape);
    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        array[index] = i;
        ++i;
    }

    array.reshape(new_shape);

    sycomore::Index index(array.dimension(), 0);
    i=0;
    for(index[1]=0; index[1]<new_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<new_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == i);
            ++i;
        }
        for(index[0]=new_shape[0]; index[0]<old_shape[0]; ++index[0])
        {
            ++i;
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeSmallerInitialized)
{
    sycomore::Shape const old_shape{7, 5};
    sycomore::Shape const new_shape{5, 3};

    sycomore::Array<float> array(old_shape);
    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        array[index] = i;
        ++i;
    }

    array.reshape(new_shape, 1000);

    sycomore::Index index(array.dimension(), 0);
    i=0;
    for(index[1]=0; index[1]<new_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<new_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == i);
            ++i;
        }
        for(index[0]=new_shape[0]; index[0]<old_shape[0]; ++index[0])
        {
            ++i;
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeMixedUninitialized)
{
    sycomore::Shape const old_shape = {7, 5};
    sycomore::Shape const new_shape = {11, 3};

    sycomore::Array<float> array(old_shape);
    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        array[index] = i;
        ++i;
    }

    array.reshape(new_shape);

    BOOST_TEST(array.shape() == new_shape);

    sycomore::Index index(array.dimension(), 0);
    i=0;
    for(index[1]=0; index[1]<new_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<old_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == i);
            ++i;
        }
    }
}

BOOST_AUTO_TEST_CASE(ReshapeMixedInitialized)
{
    sycomore::Shape const old_shape = {7, 5};
    sycomore::Shape const new_shape = {11, 3};

    sycomore::Array<float> array(old_shape);
    unsigned int i=0;
    for(auto && index: sycomore::IndexGenerator(array.shape()))
    {
        array[index] = i;
        ++i;
    }

    array.reshape(new_shape, 1000);

    BOOST_TEST(array.shape() == new_shape);

    sycomore::Index index(array.dimension(), 0);
    i=0;
    for(index[1]=0; index[1]<new_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<old_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == i);
            ++i;
        }
    }
    for(index[1]=old_shape[1]; index[1]<new_shape[1]; ++index[1])
    {
        for(index[0]=0; index[0]<new_shape[0]; ++index[0])
        {
            BOOST_TEST(array[index] == 1000);
        }
    }
}
