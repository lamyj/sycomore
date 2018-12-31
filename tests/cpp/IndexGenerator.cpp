#define BOOST_TEST_MODULE IndexGenerator
#include <boost/test/unit_test.hpp>

#include "sycomore/IndexGenerator.h"

BOOST_AUTO_TEST_CASE(WithoutOrigin)
{
    sycomore::IndexGenerator generator({3,5});
    auto && iterator = generator.begin();

    for(int y=0; y<5; ++y)
    {
        for(int x=0; x<3; ++x)
        {
            BOOST_CHECK(iterator != generator.end());

            auto && index = *iterator;
            BOOST_CHECK(index == sycomore::Index({x,y}));

            ++iterator;
        }
    }

    BOOST_CHECK(iterator == generator.end());
}

BOOST_AUTO_TEST_CASE(WithOrigin)
{
    sycomore::IndexGenerator generator({-3,-5}, {13, 11});
    auto && iterator = generator.begin();

    for(int y=-5; y<=5; ++y)
    {
        for(int x=-3; x<=9; ++x)
        {
            BOOST_CHECK(iterator != generator.end());

            auto && index = *iterator;
            BOOST_CHECK(index == sycomore::Index({x,y}));

            ++iterator;
        }
    }

    BOOST_CHECK(iterator == generator.end());
}
