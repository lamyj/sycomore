#define BOOST_TEST_MODULE Quantity
#include <boost/test/unit_test.hpp>

#include "sycomore/Quantity.h"

BOOST_AUTO_TEST_CASE(Comparison)
{
    sycomore::Quantity const q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK(q1 == q1);
    BOOST_CHECK(!(q1 == q2));
    BOOST_CHECK(!(q1 == q3));

    BOOST_CHECK(!(q1 != q1));
    BOOST_CHECK(q1 != q2);
    BOOST_CHECK(q1 != q3);

    BOOST_CHECK(q1 < q2);
    BOOST_CHECK(!(q2 <= q1));
    BOOST_CHECK(!(q1 > q2));
    BOOST_CHECK(!(q1 >= q2));
    BOOST_CHECK_THROW(q1 < q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 <= q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 > q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 >= q3, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(AdditionInPlace)
{
    sycomore::Quantity q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{5, {1,0,0,0,0,0,0}};
    q1 += q2;
    BOOST_CHECK(q1 == r);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1 += q3, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(SubtractionInPlace)
{
    sycomore::Quantity q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{-1, {1,0,0,0,0,0,0}};
    q1 -= q2;
    BOOST_CHECK(q1 == r);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1 -= q3, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(MultiplicationInPlace)
{
    sycomore::Quantity q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {0,1,0,0,0,0,0}};
    sycomore::Quantity const r{6, {1,1,0,0,0,0,0}};
    q1 *= q2;
    BOOST_CHECK(q1 == r);
}

BOOST_AUTO_TEST_CASE(ScalarMultiplicationInPlace)
{
    sycomore::Quantity q{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{6, {1,0,0,0,0,0,0}};
    q *= 3;
    BOOST_CHECK(q == r);
}

BOOST_AUTO_TEST_CASE(DivisionInPlace)
{
    sycomore::Quantity q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{4, {0,1,0,0,0,0,0}};
    sycomore::Quantity const r{0.5, {1,-1,0,0,0,0,0}};
    q1 /= q2;
    BOOST_CHECK(q1 == r);
}

BOOST_AUTO_TEST_CASE(ScalarDivisionInPlace)
{
    sycomore::Quantity q{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{0.5, {1,0,0,0,0,0,0}};
    q /= 4;
    BOOST_CHECK(q == r);
}

BOOST_AUTO_TEST_CASE(ModuloInPlace)
{
    sycomore::Quantity q1{7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{1, {1,0,0,0,0,0,0}};
    q1 %= q2;
    BOOST_CHECK(q1 == r);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1 %= q3, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ScalarModuloInPlace)
{
    sycomore::Quantity q{7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{1, {1,0,0,0,0,0,0}};
    q %= 3;
    BOOST_CHECK(q == r);
}

BOOST_AUTO_TEST_CASE(UnaryPlus)
{
    sycomore::Quantity q{2, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q == +q);
}

BOOST_AUTO_TEST_CASE(UnaryMinus)
{
    sycomore::Quantity q{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity r{-2, {1,0,0,0,0,0,0}};
    BOOST_CHECK(-q == r);
}

BOOST_AUTO_TEST_CASE(Addition)
{
    sycomore::Quantity const q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{5, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q1+q2 == r);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1+q3, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(Subtraction)
{
    sycomore::Quantity const q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{-1, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q1-q2 == r);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1-q3, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(Multiplication)
{
    sycomore::Quantity const q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {0,1,0,0,0,0,0}};
    sycomore::Quantity const r{6, {1,1,0,0,0,0,0}};
    BOOST_CHECK(q1*q2 == r);
}

BOOST_AUTO_TEST_CASE(ScalarMultiplication)
{
    sycomore::Quantity const q{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{6, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q*3 == r);
    BOOST_CHECK(3*q == r);
}

BOOST_AUTO_TEST_CASE(Division)
{
    sycomore::Quantity const q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{4, {0,1,0,0,0,0,0}};
    sycomore::Quantity const r{0.5, {1,-1,0,0,0,0,0}};
    BOOST_CHECK(q1/q2 == r);
}

BOOST_AUTO_TEST_CASE(ScalarDivision)
{
    sycomore::Quantity const q{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{0.5, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r2{1.5, {-1,0,0,0,0,0,0}};
    BOOST_CHECK(q/4 == r1);
    BOOST_CHECK(3/q == r2);
}

BOOST_AUTO_TEST_CASE(Modulo)
{
    sycomore::Quantity const q1{7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r{1, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q1%q2 == r);
}

BOOST_AUTO_TEST_CASE(ScalarModulo)
{
    sycomore::Quantity const q{7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{1, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q%3 == r1);
}

BOOST_AUTO_TEST_CASE(Power)
{
    sycomore::Quantity const q{9, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{3, {0.5,0,0,0,0,0,0}};
    BOOST_CHECK(std::pow(q, 0.5) == r1);
}
