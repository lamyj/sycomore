#define BOOST_TEST_MODULE QuantityOutOfPlaceOperators
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
        T const q1{{{1, 2}, {3, 4}}, sycomore::Length};
        T const q2{{{5, 6}, {7, 8}}, sycomore::Length};
        T const r1{{{6, 8}, {10, 12}}, sycomore::Length};
        
        sycomore::Quantity const q3{3, sycomore::Length};
        T const r2{{{4, 5}, {6, 7}}, sycomore::Length};
        
        T const q4{{{5, 6}, {7, 8}}, sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{1, sycomore::Length};
        sycomore::Quantity const q2{2, sycomore::Length};
        sycomore::Quantity const r1{3, sycomore::Length};
        
        sycomore::Quantity const q3{3, sycomore::Length};
        sycomore::Quantity const r2{4, sycomore::Length};
        
        sycomore::Quantity const q4{2, sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Addition, T, Types, Fixture<T>)
    {
        auto const t1 = this->q1 + this->q2;
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
        
        auto const t2 = this->q1 + this->q3;
        CHECK_TYPE_AND_QUANTITY(t2, T, this->r2);
        
        auto const t3 = this->q3 + this->q1;
        CHECK_TYPE_AND_QUANTITY(t3, T, this->r2);
        
        BOOST_CHECK_THROW(this->q1 + this->q4, std::runtime_error);
    }    
}

namespace Subtraction
{
    template<typename T> struct Fixture
    {
        T const q1{{{1, 2}, {3, 4}}, sycomore::Length};
        T const q2{{{8, 7}, {6, 5}}, sycomore::Length};
        T const r1{{{-7, -5}, {-3, -1}}, sycomore::Length};
        
        sycomore::Quantity const q3{3, sycomore::Length};
        T const r2{{{-2, -1}, {0, +1}}, sycomore::Length};
        T const r3{{{+2, +1}, {0, -1}}, sycomore::Length};
        
        T const q4{{{8, 7}, {6, 5}}, sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{1, sycomore::Length};
        sycomore::Quantity const q2{2, sycomore::Length};
        sycomore::Quantity const r1{-1, sycomore::Length};
        
        sycomore::Quantity const q3{3, sycomore::Length};
        sycomore::Quantity const r2{-2, sycomore::Length};
        sycomore::Quantity const r3{2, sycomore::Length};
        
        sycomore::Quantity const q4{2, sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Subtraction, T, Types, Fixture<T>)
    {
        auto const & t1 = this->q1 - this->q2;
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
        
        auto const & t2 = this->q1 - this->q3;
        CHECK_TYPE_AND_QUANTITY(t2, T, this->r2);
        
        auto const & t3 = this->q3 - this->q1;
        CHECK_TYPE_AND_QUANTITY(t3, T, this->r3);
        
        BOOST_CHECK_THROW(this->q1 - this->q4, std::runtime_error);
    }
}

namespace Multiplication
{
    template<typename T> struct Fixture
    {
        T const q1{{{2, 3}, {5, 7}}, sycomore::Length};
        T const q2{{{11, 13}, {17, 19}}, sycomore::Time};
        typename T::Container::value_type const s1{23};
        typename T::Container const s2{{29, 31}, {37, 41}};
        sycomore::Quantity const q3{3, sycomore::Time};
        
        T const r1{{{22, 39}, {85, 133}}, sycomore::Length*sycomore::Time};
        T const r2{{{46, 69}, {115, 161}}, sycomore::Length};
        T const r3{{{58, 93}, {185, 287}}, sycomore::Length};
        T const r4{{{6, 9}, {15, 21}}, sycomore::Length*sycomore::Time};
        T const r5{{{87, 93}, {111, 123}}, sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{2, sycomore::Length};
        sycomore::Quantity const q2{3, sycomore::Time};
        sycomore::Quantity::Container const s1=5;
        sycomore::Quantity::Container const s2=7;
        sycomore::Quantity const q3{3, sycomore::Time};
        
        sycomore::Quantity const r1{6, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r2{10, sycomore::Length};
        sycomore::Quantity const r3{14, sycomore::Length};
        sycomore::Quantity const r4{6, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r5{21, sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Multiplication, T, Types, Fixture<T>)
    {
        auto const t1 = this->q1 * this->q2;
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
        
        auto const t2 = this->q1 * this->s1;
        CHECK_TYPE_AND_QUANTITY(t2, T, this->r2);
        
        auto const t3 = this->q1 * this->s2;
        CHECK_TYPE_AND_QUANTITY(t3, T, this->r3);
        
        auto const t4 = this->s1 * this->q1;
        CHECK_TYPE_AND_QUANTITY(t4, T, this->r2);
        
        auto const t5 = this->s2 * this->q1;
        CHECK_TYPE_AND_QUANTITY(t5, T, this->r3);
        
        auto const t6 = this->q1 * this->q3;
        CHECK_TYPE_AND_QUANTITY(t6, T, this->r4);
        
        auto const t7 = this->q3 * this->q1;
        CHECK_TYPE_AND_QUANTITY(t7, T, this->r4);
        
        // array-of-scalar * unit
        auto const t8 = this->s2 * this->q3;
        CHECK_TYPE_AND_QUANTITY(t8, T, this->r5);
        // unit * array-of-scalar
        auto const t9 = this->q3 * this->s2;
        CHECK_TYPE_AND_QUANTITY(t9, T, this->r5);
    }
}

namespace Division
{
    template<typename T> struct Fixture
    {
        T const q1{{{14674, 27807}, {72335, 125419}}, sycomore::Length*sycomore::Time};
        T const q2{{{11, 13}, {17, 19}}, sycomore::Time};
        typename T::Container::value_type const s1{23};
        typename T::Container const s2{{29, 31}, {37, 41}};
        
        typename T::Container::value_type const s3{210};
        typename T::Container const s4{{22, 39}, {85, 133}};
        T const q3{{{2, 3}, {5, 7}}, sycomore::Time};
        sycomore::Quantity const q4{23, sycomore::Time};
        
        T const r1{{{1334, 2139}, {4255, 6601}}, sycomore::Length};
        T const r2{{{638, 1209}, {3145, 5453}}, sycomore::Length*sycomore::Time};
        T const r3{{{506, 897}, {1955, 3059}}, sycomore::Length*sycomore::Time};
        T const r4{{{105, 70}, {42, 30}}, std::pow(sycomore::Time, -1)};
        T const r5{{{11, 13}, {17, 19}}, std::pow(sycomore::Time, -1)};
        T const r6{{{638, 1209}, {3145, 5453}}, sycomore::Length};
        T const r7{
            {{29./23., 31./23.}, {37./23., 41./23.}},
            std::pow(sycomore::Time, -1)};
        T const r8{{{23./29., 23./31.}, {23./37., 23./41.}}, sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{105, sycomore::Length*sycomore::Time};
        sycomore::Quantity const q2{3, sycomore::Time};
        sycomore::Quantity::Container const s1=5;
        sycomore::Quantity::Container const s2=7;
        
        sycomore::Quantity::Container const s3{6};
        sycomore::Quantity::Container const s4{10};
        sycomore::Quantity const q3{2, sycomore::Time};
        sycomore::Quantity const q4{3, sycomore::Time};
        
        sycomore::Quantity const r1{35, sycomore::Length};
        sycomore::Quantity const r2{21, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r3{15, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r4{3, std::pow(sycomore::Time, -1)};
        sycomore::Quantity const r5{5, std::pow(sycomore::Time, -1)};
        sycomore::Quantity const r6{35, sycomore::Length};
        sycomore::Quantity const r7{7./3., std::pow(sycomore::Time, -1)};
        sycomore::Quantity const r8{3./7., sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Division, T, Types, Fixture<T>)
    {
        auto const t1 = this->q1 / this->q2;
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
        
        auto const t2 = this->q1 / this->s1;
        CHECK_TYPE_AND_QUANTITY(t2, T, this->r2);
        
        auto const t3 = this->q1 / this->s2;
        CHECK_TYPE_AND_QUANTITY(t3, T, this->r3);
        
        auto const t4 = this->s3 / this->q3;
        CHECK_TYPE_AND_QUANTITY(t4, T, this->r4);
        
        auto const t5 = this->s4 / this->q3;
        CHECK_TYPE_AND_QUANTITY(t5, T, this->r5);
        
        auto const t6 = this->q1 / this->q4;
        CHECK_TYPE_AND_QUANTITY(t6, T, this->r6);
        
        // array-of-scalar / unit
        auto const t8 = this->s2 / this->q4;
        CHECK_TYPE_AND_QUANTITY(t8, T, this->r7);
        
        // unit / array-of-scalar
        auto const t9 = this->q4 / this->s2;
        CHECK_TYPE_AND_QUANTITY(t9, T, this->r8);
    }
}

namespace Identity
{
    template<typename T> struct Fixture
    {
        T const q1{{{-1, 2}, {-3, 4}}, sycomore::Length*sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{-1, sycomore::Length*sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Identity, T, Types, Fixture<T>)
    {
        auto const t1 = +this->q1;
        CHECK_TYPE_AND_QUANTITY(t1, T, this->q1);
    }
}

namespace Opposite
{
    template<typename T> struct Fixture
    {
        T const q1{{{-1,  2}, {-3,  4}}, sycomore::Length*sycomore::Time};
        T const r1{{{ 1, -2}, { 3, -4}}, sycomore::Length*sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{-1, sycomore::Length*sycomore::Time};
        sycomore::Quantity const r1{+1, sycomore::Length*sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Opposite, T, Types, Fixture<T>)
    {
        auto const t1 = -this->q1;
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
    }
}

namespace Modulo
{
    template<typename T> struct Fixture
    {
        T const q1{{{10, 11}, {12, 13}}, sycomore::Length};
        T const q2{{{2, 3}, {8, 5}}, sycomore::Length};
        T const q3{{2, 3}, {8, 5}};
        sycomore::Quantity const q5{5, sycomore::Length};
        double const q6 = 5;
        
        typename T::Container const s1{{{2, 3}, {8, 5}}};
        
        T const r1{{{0, 2}, {4, 3}}, sycomore::Length};
        T const r2{{{0, 1}, {2, 3}}, sycomore::Length};
        
        T const q4{{{2, 3}, {8, 5}}, sycomore::Time};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{11, sycomore::Length};
        sycomore::Quantity const q2{3, sycomore::Length};
        sycomore::Quantity const q3{3};
        sycomore::Quantity const q5{5, sycomore::Length};
        double const q6 = 5;
        
        sycomore::Quantity::Container const s1{3};
        
        sycomore::Quantity const r1{2, sycomore::Length};
        sycomore::Quantity const r2{1, sycomore::Length};
        sycomore::Quantity const q4{3, sycomore::Time};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Modulo, T, Types, Fixture<T>)
    {
        auto const t1 = sycomore::fmod(this->q1, this->q2);
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
        
        auto const t2 = sycomore::fmod(this->q1, this->q3);
        CHECK_TYPE_AND_QUANTITY(t2, T, this->r1);
        
        auto const t3 = sycomore::fmod(this->q1, this->s1);
        CHECK_TYPE_AND_QUANTITY(t3, T, this->r1);
        
        auto const t4 = sycomore::fmod(this->q1, this->q5);
        CHECK_TYPE_AND_QUANTITY(t4, T, this->r2);
        
        auto const t5 = sycomore::fmod(this->q1, this->q6);
        CHECK_TYPE_AND_QUANTITY(t5, T, this->r2);
        
        BOOST_CHECK_THROW(sycomore::fmod(this->q1, this->q4), std::runtime_error);
    }
}
