#define BOOST_TEST_MODULE GridScanner
#include <boost/test/unit_test.hpp>

#include "sycomore/GridScanner.h"

sycomore::Stride compute_stride(sycomore::Shape const & shape)
{
    if(shape.size() == 0)
    {
        return sycomore::Stride();
    }

    sycomore::Stride stride(shape.size()+1);
    stride[0] = 1;
    for(unsigned int i=0; i< shape.size(); ++i)
    {
        stride[i+1] = shape[i]*stride[i];
    }

    return stride;
}

BOOST_AUTO_TEST_CASE(WithoutRegion)
{
    sycomore::Index const origin{-3, -5, -7};
    sycomore::Shape const shape{7, 11, 15};
    sycomore::Stride const stride(compute_stride(shape));

    std::vector<int> data(stride[stride.size()-1]);
    std::iota(data.begin(), data.end(), 0);

    sycomore::GridScanner const scanner(origin, shape);
    auto iterator = scanner.begin();

    for(int z=-7; z<=7; ++z)
    {
        for(int y=-5; y<=5; ++y)
        {
            for(int x=-3; x<=3; ++x)
            {
                BOOST_CHECK(iterator != scanner.end());

                auto && index = iterator->first;
                BOOST_CHECK(index == sycomore::Index({x,y, z}));

                auto && offset = iterator->second;
                BOOST_CHECK(offset == data[offset]);

                ++iterator;
            }
        }
    }
    BOOST_CHECK(iterator == scanner.end());
}

BOOST_AUTO_TEST_CASE(WithRegion3D)
{
    sycomore::Index const origin{-3, -5, -7};
    sycomore::Shape const shape{7, 11, 15};
    sycomore::Stride const stride(compute_stride(shape));

    std::vector<int> data(stride[stride.size()-1]);
    std::iota(data.begin(), data.end(), 0);

    sycomore::Index const region_origin{-1, -2, -3};
    sycomore::Shape const region_shape{5, 7, 11};

    sycomore::GridScanner const scanner(origin, shape, region_origin, region_shape);
    auto iterator = scanner.begin();

    for(int z=-3; z<=7; ++z)
    {
        for(int y=-2; y<=4; ++y)
        {
            for(int x=-1; x<=3; ++x)
            {
                BOOST_CHECK(iterator != scanner.end());

                auto && index = iterator->first;
                BOOST_CHECK(index == sycomore::Index({x,y, z}));

                auto && offset = iterator->second;
                BOOST_CHECK(offset == data[offset]);

                ++iterator;
            }
        }
    }
    BOOST_CHECK(iterator == scanner.end());

}

BOOST_AUTO_TEST_CASE(WithRegion2D)
{
    sycomore::Index const origin{-3, -5, -7};
    sycomore::Shape const shape{7, 11, 15};
    sycomore::Stride const stride(compute_stride(shape));

    std::vector<int> data(stride[stride.size()-1]);
    std::iota(data.begin(), data.end(), 0);

    sycomore::Index const region_origin{-1, -2, -3};
    sycomore::Shape const region_shape{5, 1, 11};

    sycomore::GridScanner const scanner(origin, shape, region_origin, region_shape);
    auto iterator = scanner.begin();

    for(int z=-3; z<=7; ++z)
    {
        for(int y=-2; y<=-2; ++y)
        {
            for(int x=-1; x<=3; ++x)
            {
                BOOST_CHECK(iterator != scanner.end());

                auto && index = iterator->first;
                BOOST_CHECK(index == sycomore::Index({x,y, z}));

                auto && offset = iterator->second;
                BOOST_CHECK(offset == data[offset]);

                ++iterator;
            }
        }
    }
    BOOST_CHECK(iterator == scanner.end());

}
