#define BOOST_TEST_MODULE QuantityOperators
#define BOOST_TEST_DYN_LINK

#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include "sycomore/Dimensions.h"

#include "sycomore/Quantity.h"
#include "sycomore/QuantityArray.h"
#include "sycomore/QuantityTensor.h"
#include "sycomore/QuantityTensorFixed.h"

#include "utils.h"

using Types = boost::mpl::list<
    sycomore::Quantity,
    sycomore::Matrix2x2Q, sycomore::TensorQ<2>, sycomore::ArrayQ>;

namespace Equality
{
    template<typename T> struct Fixture
    {
        T const q1{{{1, 2}, {3, 4}}, sycomore::Length};
        T const q2{{{1, 2}, {3, 4}}, sycomore::Length};
        T const q3{{{5, 6}, {7, 8}}, sycomore::Length};
        T const q4{{{1, 2}, {3, 4}}, sycomore::Time};
        T const q5{{{1, 2}, {3, 4}}, sycomore::Dimensionless};
        
        typename T::Container const scalar{{1, 2}, {3, 4}};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{1, sycomore::Length};
        sycomore::Quantity const q2{1, sycomore::Length};
        sycomore::Quantity const q3{2, sycomore::Length};
        sycomore::Quantity const q4{1, sycomore::Time};
        sycomore::Quantity const q5{2, sycomore::Dimensionless};
        
        sycomore::Quantity::Container const scalar{2};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Equal, T, Types, Fixture<T>)
    {
        BOOST_CHECK(this->q1 == this->q2);
        BOOST_CHECK(!(this->q1 == this->q3));
        BOOST_CHECK(!(this->q1 == this->q4));
        BOOST_CHECK(this->q5 == this->scalar);
        BOOST_CHECK(this->scalar == this->q5);
        BOOST_CHECK(!(this->q5 == this->q2));
    }
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Different, T, Types, Fixture<T>)
    {
        BOOST_CHECK(!(this->q1 != this->q2));
        BOOST_CHECK(this->q1 != this->q3);
        BOOST_CHECK(this->q1 != this->q4);
        BOOST_CHECK(!(this->scalar != this->q5));
        BOOST_CHECK(this->q5 != this->q2);
    }
}

BOOST_AUTO_TEST_CASE(Order)
{
    // NOTE: order is only valid for scalar quantity, not for arrays
    
    sycomore::Quantity const q1{1, sycomore::Length};
    sycomore::Quantity const q2{2, sycomore::Length};
    sycomore::Quantity const q3{1, sycomore::Time};
    
    BOOST_CHECK(q1 < q2);
    BOOST_CHECK(!(q2 < q1));
    BOOST_CHECK(!(q1 < q1));
    BOOST_CHECK_THROW(q1 < q3, std::runtime_error);
    
    BOOST_CHECK(q1 <= q2);
    BOOST_CHECK(!(q2 <= q1));
    BOOST_CHECK(q1 <= q1);
    BOOST_CHECK_THROW(q1 <= q3, std::runtime_error);
    
    BOOST_CHECK(!(q1 > q2));
    BOOST_CHECK(q2 > q1);
    BOOST_CHECK(!(q1 > q1));
    BOOST_CHECK_THROW(q1 > q3, std::runtime_error);
    
    BOOST_CHECK(!(q1 >= q2));
    BOOST_CHECK(q2 >= q1);
    BOOST_CHECK(q1 >= q1);
    BOOST_CHECK_THROW(q1 >= q3, std::runtime_error);
}
