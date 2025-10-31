#define BOOST_TEST_MODULE QuantityInPlaceOperators
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

namespace Addition
{
    template<typename T> struct Fixture
    {
        T q1{{{1, 2}, {3, 4}}, sycomore::Length};
        T const q2{{{5, 6}, {7, 8}}, sycomore::Length};
        T const r1{{{6, 8}, {10, 12}}, sycomore::Length};
        
        sycomore::Quantity const q4{3, sycomore::Length};
        T const r2{{{9, 11}, {13, 15}}, sycomore::Length};
        
        T const q3{{{5, 6}, {7, 8}}, sycomore::Time};
        
        T q5{{{1, 2}, {3, 4}}, sycomore::Dimensionless};
        typename T::Container const scalar{{5, 6}, {7, 8}};
        T const r3{{{6, 8}, {10, 12}}, sycomore::Dimensionless};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity q1{1, sycomore::Length};
        sycomore::Quantity const q2{2, sycomore::Length};
        sycomore::Quantity const r1{3, sycomore::Length};
        
        sycomore::Quantity const q4{3, sycomore::Length};
        sycomore::Quantity const r2{6, sycomore::Length};
        
        sycomore::Quantity const q3{2, sycomore::Time};
        
        sycomore::Quantity q5{1, sycomore::Dimensionless};
        sycomore::Quantity::Container const scalar{2};
        sycomore::Quantity const r3{3, sycomore::Dimensionless};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Addition, T, Types, Fixture<T>)
    {
        auto const & result1 = (this->q1 += this->q2);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, result1, this->r1);
        
        BOOST_CHECK_THROW(this->q1 += this->q3, std::runtime_error);
        
        auto const & result2 = (this->q1 += this->q4);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, result2, this->r2);
        
        auto const & result3 = (this->q5 += this->scalar);
        CHECK_IDENTITY_AND_QUANTITY(this->q5, result3, this->r3);
    }
}

namespace Subtraction
{
    template<typename T> struct Fixture
    {
        T q1{{{1, 2}, {3, 4}}, sycomore::Length};
        T const q2{{{8, 7}, {6, 5}}, sycomore::Length};
        T const r1{{{-7, -5}, {-3, -1}}, sycomore::Length};
        
        T const q3{{{8, 7}, {6, 5}}, sycomore::Time};
        
        sycomore::Quantity const q4{3, sycomore::Length};
        T const r2{{{-10, -8}, {-6, -4}}, sycomore::Length};
        
        T q5{{{1, 2}, {3, 4}}, sycomore::Dimensionless};
        typename T::Container const scalar{{8, 7}, {6, 5}};
        T const r3{{{-7, -5}, {-3, -1}}, sycomore::Dimensionless};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity q1{1, sycomore::Length};
        sycomore::Quantity const q2{2, sycomore::Length};
        sycomore::Quantity const r1{-1, sycomore::Length};
        
        sycomore::Quantity const q3{2, sycomore::Time};
        
        sycomore::Quantity const q4{3, sycomore::Length};
        sycomore::Quantity const r2{-4, sycomore::Length};
        
        sycomore::Quantity q5{1, sycomore::Dimensionless};
        sycomore::Quantity::Container const scalar{2};
        sycomore::Quantity const r3{-1, sycomore::Dimensionless};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Subtraction, T, Types, Fixture<T>)
    {
        auto const & result = (this->q1 -= this->q2);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, result, this->r1);
        
        BOOST_CHECK_THROW(this->q1 -= this->q3, std::runtime_error);
        
        auto const & result2 = (this->q1 -= this->q4);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, result2, this->r2);
        
        auto const & result3 = (this->q5 -= this->scalar);
        CHECK_IDENTITY_AND_QUANTITY(this->q5, result3, this->r3);
    }
}

namespace Multiplication
{
    template<typename T> struct Fixture
    {
        T q1{{{2, 3}, {5, 7}}, sycomore::Length};
        T const q2{{{11, 13}, {17, 19}}, sycomore::Time};
        typename T::Container::value_type const s1{23};
        typename T::Container const s2{{29, 31}, {37, 41}};
        sycomore::Quantity const q3{3, sycomore::Length};
        
        T const r1{{{22, 39}, {85, 133}}, sycomore::Length*sycomore::Time};
        T const r2{{{506, 897}, {1955, 3059}}, sycomore::Length*sycomore::Time};
        T const r3{
            {{14674, 27807}, {72335, 125419}}, sycomore::Length*sycomore::Time};
        T const r4{
            {{44022, 83421}, {217005, 376257}},
            std::pow(sycomore::Length, 2)*sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity q1{2, sycomore::Length};
        sycomore::Quantity const q2{3, sycomore::Time};
        sycomore::Quantity::Container const s1=5;
        sycomore::Quantity::Container const s2=7;
        sycomore::Quantity const q3{3, sycomore::Length};
        
        sycomore::Quantity const r1{6, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r2{30, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r3{210, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r4{630, std::pow(sycomore::Length, 2)*sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Multiplication, T, Types, Fixture<T>)
    {
        auto const & t1 = (this->q1 *= this->q2);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t1, this->r1);
        
        auto const & t2 = (this->q1 *= this->s1);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t2, this->r2);
        
        auto const & t3 = (this->q1 *= this->s2);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t3, this->r3);
    }
}

namespace Division
{
    template<typename T> struct Fixture
    {
        T q1{{{14674, 27807}, {72335, 125419}}, sycomore::Length};
        T const q2{{{29, 31}, {37, 41}}, sycomore::Time};
        typename T::Container::value_type const s1{23};
        typename T::Container const s2{{11, 13}, {17, 19}};
        
        T const r1{{{506, 897}, {1955, 3059}}, sycomore::Length/sycomore::Time};
        T const r2{{{22, 39}, {85, 133}}, sycomore::Length/sycomore::Time};
        T const r3{{{2, 3}, {5, 7}}, sycomore::Length/sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity q1{210, sycomore::Length};
        sycomore::Quantity const q2{7, sycomore::Time};
        sycomore::Quantity::Container const s1=5;
        sycomore::Quantity::Container const s2=2;
        
        sycomore::Quantity const r1{30, sycomore::Length/sycomore::Time};
        sycomore::Quantity const r2{6, sycomore::Length/sycomore::Time};
        sycomore::Quantity const r3{3, sycomore::Length/sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Division, T, Types, Fixture<T>)
    {
        auto const & t1 = (this->q1 /= this->q2);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t1, this->r1);
        
        auto const & t2 = (this->q1 /= this->s1);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t2, this->r2);
        
        auto const & t3 = (this->q1 /= this->s2);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t3, this->r3);
    }
}

namespace Modulo
{
    template<typename T> struct Fixture
    {
        T q1{{{10, 11}, {12, 13}}, sycomore::Length};
        T const q2{{{2, 3}, {8, 5}}, sycomore::Length};
        T const q3{{{2, 3}, {8, 5}}, sycomore::Dimensionless};
        
        typename T::Container const s1{{{2, 3}, {8, 5}}};
        
        T const r1{{{0, 2}, {4, 3}}, sycomore::Length};
        
        T const q4{{{2, 3}, {8, 5}}, sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity q1{11, sycomore::Length};
        sycomore::Quantity const q2{3, sycomore::Length};
        sycomore::Quantity const q3{3};
        
        sycomore::Quantity::Container const s1{3};
        
        sycomore::Quantity const r1{2, sycomore::Length};
        
        sycomore::Quantity const q4{3, sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(ModuloQuantity, T, Types, Fixture<T>)
    {
        auto const & t = (this->q1 %= this->q2);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t, this->r1);
        
        BOOST_CHECK_THROW(this->q1 %= this->q4, std::runtime_error);
    }
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(ModuloScalarQuantity, T, Types, Fixture<T>)
    {
        auto const & t = (this->q1 %= this->q3);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t, this->r1);
    }
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(ModuloScalarArray, T, Types, Fixture<T>)
    {
        auto const & t = (this->q1 %= this->s1);
        CHECK_IDENTITY_AND_QUANTITY(this->q1, t, this->r1);
    }
}
