#define BOOST_TEST_MODULE sycomore::QuantityAccess
#define BOOST_TEST_DYN_LINK

#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include "sycomore/Dimensions.h"

#include "sycomore/Quantity.h"
#include "sycomore/QuantityArray.h"
#include "sycomore/QuantityConstIterator.h"
#include "sycomore/QuantityIterator.h"
#include "sycomore/QuantityTensor.h"
#include "sycomore/QuantityTensorFixed.h"

#include "utils.h"

using Types = boost::mpl::list<
    sycomore::Matrix2x2Q, sycomore::TensorQ<2>, sycomore::ArrayQ>;

template<typename T> struct Fixture
{
    T const q1{{{1, 2}, {3, 4}}, sycomore::Length};
    T q2{{{1, 2}, {3, 4}}, sycomore::Length};
};

BOOST_FIXTURE_TEST_CASE_TEMPLATE(ConstIterator, T, Types, Fixture<T>)
{
    sycomore::QuantityConstIterator<T> iterator(this->q1);
    sycomore::QuantityConstIterator<T> other(this->q1);
    sycomore::QuantityConstIterator<T> const end(this->q1, true);
    
    BOOST_CHECK(iterator == other);
    BOOST_CHECK(!(iterator == end));
    
    BOOST_CHECK(!(iterator != other));
    BOOST_CHECK(iterator != end);
    
    CHECK_QUANTITY(*iterator, sycomore::Quantity(1, sycomore::Length));
    
    auto const r1 = *(++iterator);
    CHECK_QUANTITY(r1, sycomore::Quantity(2, sycomore::Length));
    
    auto const it2 = iterator++;
    CHECK_QUANTITY(*it2, sycomore::Quantity(2, sycomore::Length));
    CHECK_QUANTITY(*iterator, sycomore::Quantity(3, sycomore::Length));
    
    auto const r3 = *(--iterator);
    CHECK_QUANTITY(r3, sycomore::Quantity(2, sycomore::Length));
    
    auto const it3 = iterator--;
    CHECK_QUANTITY(*it3, sycomore::Quantity(2, sycomore::Length));
    CHECK_QUANTITY(*iterator, sycomore::Quantity(1, sycomore::Length));
    
    iterator = sycomore::QuantityConstIterator<T>(this->q1);
    ++iterator;
    ++iterator;
    ++iterator;
    ++iterator;
    BOOST_CHECK(iterator == end);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(Iterator, T, Types, Fixture<T>)
{
    sycomore::QuantityIterator<T> iterator(this->q2);
    sycomore::QuantityIterator<T> other(this->q2);
    sycomore::QuantityIterator<T> const end(this->q2, true);
    
    BOOST_CHECK(iterator == other);
    BOOST_CHECK(!(iterator == end));
    
    BOOST_CHECK(!(iterator != other));
    BOOST_CHECK(iterator != end);
    
    CHECK_QUANTITY(*iterator, sycomore::Quantity(1, sycomore::Length));
    
    auto const r1 = *(++iterator);
    CHECK_QUANTITY(r1, sycomore::Quantity(2, sycomore::Length));
    
    auto const it2 = iterator++;
    CHECK_QUANTITY(*it2, sycomore::Quantity(2, sycomore::Length));
    CHECK_QUANTITY(*iterator, sycomore::Quantity(3, sycomore::Length));
    
    auto const r3 = *(--iterator);
    CHECK_QUANTITY(r3, sycomore::Quantity(2, sycomore::Length));
    
    auto const it3 = iterator--;
    CHECK_QUANTITY(*it3, sycomore::Quantity(2, sycomore::Length));
    CHECK_QUANTITY(*iterator, sycomore::Quantity(1, sycomore::Length));
    
    (*iterator).magnitude = 42;
    CHECK_QUANTITY(*iterator, sycomore::Quantity(42, sycomore::Length));
    
    (*iterator) = sycomore::Quantity(43, sycomore::Length);
    CHECK_QUANTITY(*iterator, sycomore::Quantity(43, sycomore::Length));
    
    BOOST_CHECK_THROW(
        (*iterator) = sycomore::Quantity(43, sycomore::Time),
        std::runtime_error);
    
    iterator = sycomore::QuantityIterator<T>(this->q2);
    ++iterator;
    ++iterator;
    ++iterator;
    ++iterator;
    BOOST_CHECK(iterator == end);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(BeginEnd, T, Types, Fixture<T>)
{
    CHECK_QUANTITY(*this->q1.begin(), sycomore::Quantity(1, sycomore::Length));
    CHECK_QUANTITY(*(--this->q1.end()), sycomore::Quantity(4, sycomore::Length));
    
    CHECK_QUANTITY(*this->q1.cbegin(), sycomore::Quantity(1, sycomore::Length));
    CHECK_QUANTITY(*(--this->q1.cend()), sycomore::Quantity(4, sycomore::Length));
    
    CHECK_QUANTITY(*this->q2.begin(), sycomore::Quantity(1, sycomore::Length));
    CHECK_QUANTITY(*(--this->q2.end()), sycomore::Quantity(4, sycomore::Length));
    
    *this->q2.begin() = sycomore::Quantity(42, sycomore::Length);
    CHECK_QUANTITY(*this->q2.begin(), sycomore::Quantity(42, sycomore::Length));
    
    int magnitude = 1;
    for(auto && q: this->q1)
    {
        CHECK_QUANTITY(q, sycomore::Quantity(magnitude, sycomore::Length));
        ++magnitude;
    }
}
