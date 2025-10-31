#define BOOST_TEST_MODULE QuantityMisc
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

BOOST_AUTO_TEST_CASE_TEMPLATE(IsQuantity, T, Types)
{
    BOOST_CHECK(sycomore::is_quantity<T>::value);
    BOOST_CHECK(!sycomore::is_quantity<typename T::Container>::value);
}

namespace Conversion
{
    struct FixtureBase
    {
        sycomore::Quantity const mm{1e-3, sycomore::Length};
        sycomore::Quantity const s{1, sycomore::Time};
    };
    template<typename T> struct Fixture: public FixtureBase
    {
        T const q1{{{1, 2}, {3, 4}}, sycomore::Length};
        typename T::Container const r1{{1000, 2000}, {3000, 4000}};
    };
    
    template<> struct Fixture<sycomore::Quantity>: public FixtureBase
    {
        sycomore::Quantity const q1{2, sycomore::Length};
        sycomore::Quantity::Container const r1 = 2000;
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Addition, T, Types, Fixture<T>)
    {
        auto const t = this->q1.convert_to(this->mm);
        BOOST_CHECK(typeid(t) == typeid(typename T::Container));
        BOOST_CHECK(t == this->r1);
        
        BOOST_CHECK_THROW(this->q1.convert_to(this->s), std::runtime_error);
    }
}

namespace Abs
{
    template<typename T> struct Fixture
    {
        T const q1{{{-1, 2}, {-3,  4}}, sycomore::Length};
        sycomore::QuantityConstReference const q2{q1(0, 0)};
        T const r1{{{ 1, 2}, { 3, 4}}, sycomore::Length};
        sycomore::Quantity const r2{r1(0, 0).magnitude, r1.dimensions};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{-1, sycomore::Length};
        sycomore::Quantity const q2{q1};
        sycomore::Quantity const r1{+1, sycomore::Length};
        sycomore::Quantity const r2{r1};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Abs, T, Types, Fixture<T>)
    {
        auto const t1 = std::abs(this->q1);
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
        
        auto const t2 = std::abs(this->q2);
        CHECK_TYPE_AND_QUANTITY(t2, sycomore::Quantity, this->r2);
    }
}

namespace Pow
{
    template<typename T> struct Fixture
    {
        T const q1{{{-1, 2}, { -3,  4}}, sycomore::Length};
        T const r1{{{-1, 8}, {-27, 64}}, std::pow(sycomore::Length, 3)};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{-2, sycomore::Length};
        sycomore::Quantity const r1{-8, std::pow(sycomore::Length, 3)};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Pow, T, Types, Fixture<T>)
    {
        auto const t1 = std::pow(this->q1, 3);
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
    }
}

namespace RoundAndFriends
{
    template<typename T> struct Fixture
    {
        T const q1{{{ 9.2, - 9.7}, { -9,  0}}, std::pow(sycomore::Length, 1.5)};
        T const r1{{{ 9  , -10  }, { -9,  0}}, std::pow(sycomore::Length, 1.5)};
        T const r2{{{ 9  , - 9  }, { -9,  0}}, std::pow(sycomore::Length, 1.5)};
        T const r3{{{ 9  , -10  }, { -9,  0}}, std::pow(sycomore::Length, 1.5)};
        T const r4{{{10  , - 9  }, { -9,  0}}, std::pow(sycomore::Length, 1.5)};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q1{- 9.7, std::pow(sycomore::Length, 1.5)};
        sycomore::Quantity const r1{-10  , std::pow(sycomore::Length, 1.5)};
        sycomore::Quantity const r2{- 9  , std::pow(sycomore::Length, 1.5)};
        sycomore::Quantity const r3{-10  , std::pow(sycomore::Length, 1.5)};
        sycomore::Quantity const r4{- 9  , std::pow(sycomore::Length, 1.5)};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(RoundAndFriends, T, Types, Fixture<T>)
    {
        auto const t1 = std::round(this->q1);
        CHECK_TYPE_AND_QUANTITY(t1, T, this->r1);
        
        auto const t2 = std::trunc(this->q1);
        CHECK_TYPE_AND_QUANTITY(t2, T, this->r2);
        
        auto const t3 = std::floor(this->q1);
        CHECK_TYPE_AND_QUANTITY(t3, T, this->r3);
        
        auto const t4 = std::ceil(this->q1);
        CHECK_TYPE_AND_QUANTITY(t4, T, this->r4);
    }
}

namespace Hash
{
    template<typename T> struct Fixture
    {
        T const q{{{ 9.2, - 9.7}, { -9,  0}}, std::pow(sycomore::Length, 1.5)};
    };
    
    template<> struct Fixture<sycomore::Quantity>
    {
        sycomore::Quantity const q{- 9.7, std::pow(sycomore::Length, 1.5)};
    };
    
    BOOST_FIXTURE_TEST_CASE_TEMPLATE(Hash, T, Types, Fixture<T>)
    {
        BOOST_CHECK(std::hash<T>{}(this->q) != 0);
    }
}
