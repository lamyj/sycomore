#define BOOST_TEST_MODULE isochromat_Operator
#include <boost/test/unit_test.hpp>

#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#include "sycomore/isochromat/Operator.h"

#include <iostream>
#include <xtensor/xio.hpp>

BOOST_AUTO_TEST_CASE(PostMultiply)
{
    sycomore::isochromat::Operator const left({
        {
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9, 10, 11, 12},
            {13, 14, 15, 16}
        },
        {
            {13, 9, 5, 1},
            {14, 10, 6, 2},
            {15, 11, 7, 3},
            {16, 12, 8, 4}
        }
    });
    
    sycomore::isochromat::Operator right({
        {
            {16, 15, 14, 13},
            {12, 11, 10, 9},
            {8, 7, 6, 5},
            {4, 3, 2, 1}
        },
        {
            {4, 3, 2, 1},
            {8, 7, 6, 5},
            {12, 11, 10, 9},
            {16, 15, 14, 13}
        }
    });
    
    right.preMultiply(left);
    
    sycomore::isochromat::Operator::Array expected{
        {{80, 70, 60, 50},
         {240, 214, 188, 162},
         {400, 358, 316, 274},
         {560, 502, 444, 386}},
        {{200, 172, 144, 116},
         {240, 208, 176, 144},
         {280, 244, 208, 172},
         {320, 280, 240, 200}}
    };
    BOOST_TEST(xt::allclose(right.array(), expected));
}

BOOST_AUTO_TEST_CASE(PreMultiply)
{
    sycomore::isochromat::Operator left({
        {
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9, 10, 11, 12},
            {13, 14, 15, 16}
        },
        {
            {13, 9, 5, 1},
            {14, 10, 6, 2},
            {15, 11, 7, 3},
            {16, 12, 8, 4}
        }
    });
    
    sycomore::isochromat::Operator const right({
        {
            {16, 15, 14, 13},
            {12, 11, 10, 9},
            {8, 7, 6, 5},
            {4, 3, 2, 1}
        },
        {
            {4, 3, 2, 1},
            {8, 7, 6, 5},
            {12, 11, 10, 9},
            {16, 15, 14, 13}
        }
    });
    
    left *= right;
    
    sycomore::isochromat::Operator::Array expected{
        {{80, 70, 60, 50},
         {240, 214, 188, 162},
         {400, 358, 316, 274},
         {560, 502, 444, 386}},
        {{200, 172, 144, 116},
         {240, 208, 176, 144},
         {280, 244, 208, 172},
         {320, 280, 240, 200}}
    };
    BOOST_TEST(xt::allclose(left.array(), expected));
}

BOOST_AUTO_TEST_CASE(Multiply)
{
    sycomore::isochromat::Operator left({
        {
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9, 10, 11, 12},
            {13, 14, 15, 16}
        },
        {
            {13, 9, 5, 1},
            {14, 10, 6, 2},
            {15, 11, 7, 3},
            {16, 12, 8, 4}
        }
    });
    
    sycomore::isochromat::Operator const right({
        {
            {16, 15, 14, 13},
            {12, 11, 10, 9},
            {8, 7, 6, 5},
            {4, 3, 2, 1}
        },
        {
            {4, 3, 2, 1},
            {8, 7, 6, 5},
            {12, 11, 10, 9},
            {16, 15, 14, 13}
        }
    });
    
    auto const combined = left * right;
    
    sycomore::isochromat::Operator::Array expected{
        {{80, 70, 60, 50},
         {240, 214, 188, 162},
         {400, 358, 316, 274},
         {560, 502, 444, 386}},
        {{200, 172, 144, 116},
         {240, 208, 176, 144},
         {280, 244, 208, 172},
         {320, 280, 240, 200}}
    };
    BOOST_TEST(xt::allclose(combined.array(), expected));
}
