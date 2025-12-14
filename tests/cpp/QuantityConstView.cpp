#define BOOST_TEST_MODULE QuantityConstView
#define BOOST_TEST_DYN_LINK

#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include "sycomore/Dimensions.h"

#include "sycomore/QuantityArray.h"
#include "sycomore/QuantityConstView.h"
#include "sycomore/QuantityTensor.h"

#include "utils.h"

using Types = boost::mpl::list<sycomore::TensorQ<3>, sycomore::ArrayQ>;

namespace OneDimension
{
    template<typename T> struct Fixture
    {
        T const q{
            {
                { { 1,  2,  3,  4}, { 5,  6,  7,  8} },
                { { 9, 10, 11, 12}, {13, 14, 15, 16} },
                { {17, 18, 19, 20}, {21, 22, 23, 24} },
            },
            sycomore::Length};
        
        sycomore::TensorQ<2> const t{
            { { 9, 10, 11, 12}, {13, 14, 15, 16} }, sycomore::Length};
        sycomore::ArrayQ const a{
            { { 9, 10, 11, 12}, {13, 14, 15, 16} }, sycomore::Length};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(OneDimension, T, Types, Fixture<T>)
    {
        auto const v = sycomore::view(this->q, 1);
        
        BOOST_CHECK(v.size() == 8);
        
        auto const v_shape = v.shape();
        std::array<int, 2> expected_shape{2, 4};
        BOOST_CHECK(v_shape.size() == expected_shape.size());
        BOOST_CHECK(std::equal(v_shape.begin(), v_shape.end(), expected_shape.begin()));
        for(std::size_t i=0; i!=v_shape.size(); ++i)
        {
            BOOST_CHECK(v.shape(i) == expected_shape[i]);
        }
        
        CHECK_QUANTITY(v(1, 2), sycomore::Quantity(15, sycomore::Length));
        CHECK_QUANTITY(v.at(1, 2), sycomore::Quantity(15, sycomore::Length));
        CHECK_QUANTITY(v.unchecked(1, 2), sycomore::Quantity(15, sycomore::Length));
        
        std::array<int, 2> index{1, 2};
        CHECK_QUANTITY(v[index], sycomore::Quantity(15, sycomore::Length));
        CHECK_QUANTITY((v[{1, 2}]), sycomore::Quantity(15, sycomore::Length));
        
        CHECK_QUANTITY(sycomore::TensorQ<2>(v), this->t);
        CHECK_QUANTITY(sycomore::ArrayQ(v), this->a);
    }
}

namespace All
{
    template<typename T> struct Fixture
    {
        T const q{
            {
                { { 1,  2,  3,  4}, { 5,  6,  7,  8} },
                { { 9, 10, 11, 12}, {13, 14, 15, 16} },
                { {17, 18, 19, 20}, {21, 22, 23, 24} },
            },
            sycomore::Length};
        
        sycomore::TensorQ<2> const t{
            { { 5,  6,  7,  8}, {13, 14, 15, 16}, {21, 22, 23, 24} },
            sycomore::Length};
        sycomore::ArrayQ const a{
            { { 5,  6,  7,  8}, {13, 14, 15, 16}, {21, 22, 23, 24} },
            sycomore::Length};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(All, T, Types, Fixture<T>)
    {
        auto const v = sycomore::view(this->q, xt::all(), 1);
        CHECK_QUANTITY(sycomore::TensorQ<2>(v), this->t);
        CHECK_QUANTITY(sycomore::ArrayQ(v), this->a);
    }
}

namespace Range
{
    template<typename T> struct Fixture
    {
        T const q{
            {
                { { 1,  2,  3,  4}, { 5,  6,  7,  8} },
                { { 9, 10, 11, 12}, {13, 14, 15, 16} },
                { {17, 18, 19, 20}, {21, 22, 23, 24} },
            },
            sycomore::Length};
        
        sycomore::TensorQ<3> const t{
            { 
                { { 2,  3}, { 6,  7} },
                { {10, 11}, {14, 15} },
                { {18, 19}, {22, 23} }
            },
            sycomore::Length};
        sycomore::ArrayQ const a{
            { 
                { { 2,  3}, { 6,  7} },
                { {10, 11}, {14, 15} },
                { {18, 19}, {22, 23} }
            },
            sycomore::Length};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Range, T, Types, Fixture<T>)
    {
        auto const v = sycomore::view(this->q, xt::all(), xt::all(), xt::range(1, 3));
        CHECK_QUANTITY(sycomore::TensorQ<3>(v), this->t);
        CHECK_QUANTITY(sycomore::ArrayQ(v), this->a);
    }
}

namespace Scalar
{
    template<typename T> struct Fixture
    {
        T const q{
            {
                { { 1,  2,  3,  4}, { 5,  6,  7,  8} },
                { { 9, 10, 11, 12}, {13, 14, 15, 16} },
                { {17, 18, 19, 20}, {21, 22, 23, 24} },
            },
            sycomore::Length};
        
        sycomore::Quantity const r{18, sycomore::Length};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Scalar, T, Types, Fixture<T>)
    {
        auto const v = sycomore::view(this->q, 2, 0, 1);
        CHECK_QUANTITY(sycomore::Quantity(v), this->r);
    }
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(ScalarError, T, Types, Fixture<T>)
    {
        BOOST_CHECK_THROW(sycomore::Quantity(sycomore::view(this->q, 2, 0)), std::runtime_error);
    }
}

namespace Strided
{
    template<typename T> struct Fixture
    {
        T const q{
            {
                { { 1,  2,  3,  4}, { 5,  6,  7,  8} },
                { { 9, 10, 11, 12}, {13, 14, 15, 16} },
                { {17, 18, 19, 20}, {21, 22, 23, 24} },
            },
            sycomore::Length};
        
        sycomore::TensorQ<3> const t{
            { 
                { { 2,  3}, { 6,  7} },
                { {10, 11}, {14, 15} },
                { {18, 19}, {22, 23} }
            },
            sycomore::Length};
        sycomore::ArrayQ const a{
            { 
                { { 2,  3}, { 6,  7} },
                { {10, 11}, {14, 15} },
                { {18, 19}, {22, 23} }
            },
            sycomore::Length};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Strided, T, Types, Fixture<T>)
    {
        auto const v = sycomore::strided_view(
            this->q, {xt::all(), xt::all(), xt::range(1, 3)});
        CHECK_QUANTITY(sycomore::TensorQ<3>(v), this->t);
        CHECK_QUANTITY(sycomore::ArrayQ(v), this->a);
    }
}
