#define BOOST_TEST_MODULE sycomore::QuantityAccess
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
    sycomore::Matrix2x2Q, sycomore::TensorQ<2>, sycomore::ArrayQ>;

template<typename T> struct Fixture
{
    T const q1{{{1, 2}, {3, 4}}, sycomore::Length};
    T q2{{{1, 2}, {3, 4}}, sycomore::Length};
};

BOOST_FIXTURE_TEST_CASE_TEMPLATE(CallOperator, T, Types, Fixture<T>)
{
    using sycomore::Quantity;
    using sycomore::Length;
    using sycomore::Time;
    
    CHECK_IDENTITY(this->q1(1, 1).magnitude, this->q1.magnitude(1, 1));
    CHECK_QUANTITY(
        this->q1(1, 1), Quantity(this->q1.magnitude(1, 1), this->q1.dimensions));
    
    this->q2(1, 1) = {42., Length};
    CHECK_QUANTITY(this->q2(1, 1), Quantity(42, this->q2.dimensions));
    
    BOOST_CHECK_THROW(
        (this->q2(1, 1) = Quantity{42., Time}), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(At, T, Types, Fixture<T>)
{
    using sycomore::Quantity;
    using sycomore::Length;
    using sycomore::Time;
    
    CHECK_IDENTITY(this->q1.at(1, 1).magnitude, this->q1.magnitude.at(1, 1));
    CHECK_QUANTITY(
        this->q1.at(1, 1),
        Quantity(this->q1.magnitude.at(1, 1), this->q1.dimensions));
    
    this->q2.at(1, 1) = {42., Length};
    CHECK_IDENTITY(this->q1.at(1, 1).magnitude, this->q1.magnitude.at(1, 1));
    CHECK_QUANTITY(this->q2.at(1, 1), Quantity(42, this->q2.dimensions));
    
    BOOST_CHECK_THROW(
        (this->q2.at(1, 1) = Quantity{42., Time}), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(Unchecked, T, Types, Fixture<T>)
{
    using sycomore::Quantity;
    using sycomore::Length;
    using sycomore::Time;
    
    CHECK_IDENTITY(
        this->q1.unchecked(1, 1).magnitude, this->q1.magnitude.unchecked(1, 1));
    CHECK_QUANTITY(
        this->q1.unchecked(1, 1),
        Quantity(this->q1.magnitude.unchecked(1, 1), this->q1.dimensions));
    
    this->q2.unchecked(1, 1) = {42., Length};
    CHECK_IDENTITY(
        this->q2.unchecked(1, 1).magnitude, this->q2.magnitude.unchecked(1, 1));
    CHECK_QUANTITY(this->q2.unchecked(1, 1), Quantity(42, this->q2.dimensions));
    
    BOOST_CHECK_THROW(
        (this->q2.unchecked(1, 1) = Quantity{42., Time}), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(SubscriptOperator, T, Types, Fixture<T>)
{
    using sycomore::Quantity;
    using sycomore::Length;
    using sycomore::Time;
    
    CHECK_IDENTITY((this->q1[{1, 1}].magnitude), (this->q1.magnitude[{1, 1}]));
    CHECK_QUANTITY(
        (this->q1[{1, 1}]),
        (Quantity(this->q1.magnitude[{1, 1}], this->q1.dimensions)));
    
    this->q2[{1, 1}] = {42., Length};
    CHECK_IDENTITY((this->q2[{1, 1}].magnitude), (this->q2.magnitude[{1, 1}]));
    CHECK_QUANTITY((this->q2[{1, 1}]), Quantity(42, this->q2.dimensions));
    
    BOOST_CHECK_THROW(
        (this->q2[{1, 1}] = Quantity{42., Time}), std::runtime_error);
}

template<typename T> struct QuantityReferenceFixture: public Fixture<T>
{
    template<typename T2>
    static void test(T2 && q)
    {
        using sycomore::Quantity;
        using sycomore::Length;
        using sycomore::Time;
        
        // Addition
        CHECK_QUANTITY(q+Quantity(3, Length), Quantity(7, Length));
        CHECK_QUANTITY(Quantity(3, Length)+q, Quantity(7, Length));
        // Subtraction
        CHECK_QUANTITY(q-Quantity(1, Length), Quantity(3, Length));
        CHECK_QUANTITY(Quantity(1, Length)-q, Quantity(-3, Length));
        // Multiplication
        CHECK_QUANTITY(q*2, Quantity(8, Length));
        CHECK_QUANTITY(2*q, Quantity(8, Length));
        CHECK_QUANTITY(q*Quantity(3, Time), Quantity(12, Length*Time));
        CHECK_QUANTITY(Quantity(3, Time)*q, Quantity(12, Length*Time));
        // Division
        CHECK_QUANTITY(q/.5, Quantity(8, Length));
        CHECK_QUANTITY(8/q, Quantity(2, std::pow(Length, -1)));
        CHECK_QUANTITY(q/Quantity(8, Time), Quantity(0.5, Length/Time));
        CHECK_QUANTITY(Quantity(8, Time)/q, Quantity(2, Time/Length));
        // Modulo
        CHECK_QUANTITY(sycomore::fmod(q, 3.), Quantity(1, Length));
        CHECK_QUANTITY(sycomore::fmod(q, Quantity(3, Length)), Quantity(1, Length));
        CHECK_QUANTITY(sycomore::fmod(Quantity(5, Length), q), Quantity(1, Length));
        // Identity
        CHECK_QUANTITY(+q, Quantity(+4, Length));
        // Opposite
        CHECK_QUANTITY(-q, Quantity(-4, Length));
    }

};

BOOST_FIXTURE_TEST_CASE_TEMPLATE(QuantityConstReference, T, Types, QuantityReferenceFixture<T>)
{
    this->test(this->q1(1, 1));
    auto const q2 = 3 * this->q1(1, 1);
    CHECK_QUANTITY(q2, sycomore::Quantity(12, sycomore::Length));
    auto const q3 = this->q1(1, 1) + sycomore::Quantity(3, sycomore::Length);
    CHECK_QUANTITY(q3, sycomore::Quantity(7, sycomore::Length));
    BOOST_CHECK_THROW(
        this->q1(1, 1) + sycomore::Quantity(3, sycomore::Time),
        std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(QuantityReference, T, Types, QuantityReferenceFixture<T>)
{
    this->test(this->q2(1, 1));
    this->q2(1, 1) *= 3;
    CHECK_QUANTITY(this->q2(1, 1), sycomore::Quantity(12, sycomore::Length));
    this->q2(1, 1) += sycomore::Quantity(3, sycomore::Length);
    CHECK_QUANTITY(this->q2(1, 1), sycomore::Quantity(15, sycomore::Length));
    BOOST_CHECK_THROW(
        this->q2(1, 1) += sycomore::Quantity(3, sycomore::Time),
        std::runtime_error);
}
