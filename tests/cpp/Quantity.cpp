#define BOOST_TEST_MODULE Quantity
#include <boost/test/unit_test.hpp>

#include "sycomore/Quantity.h"

BOOST_AUTO_TEST_CASE(Comparison)
{
    sycomore::Quantity const q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    sycomore::Quantity const scalar{2, {0,0,0,0,0,0,0}};
    
    BOOST_CHECK(q1 == q1);
    BOOST_CHECK(!(q1 == q2));
    BOOST_CHECK(!(q1 == q3));
    BOOST_CHECK(scalar == 2);
    BOOST_CHECK(2 == scalar);
    BOOST_CHECK(!(scalar == 3));
    BOOST_CHECK(!(3 == scalar));

    BOOST_CHECK(!(q1 != q1));
    BOOST_CHECK(q1 != q2);
    BOOST_CHECK(q1 != q3);
    BOOST_CHECK(!(scalar != 2));
    BOOST_CHECK(!(2 != scalar));
    BOOST_CHECK(scalar != 3);
    BOOST_CHECK(3 !=scalar);

    BOOST_CHECK(q1 < q2);
    BOOST_CHECK(scalar < 3);
    BOOST_CHECK(1 < scalar);
    BOOST_CHECK(!(q2 <= q1));
    BOOST_CHECK(scalar <= 2);
    BOOST_CHECK(2 <= scalar);
    BOOST_CHECK(!(q1 > q2));
    BOOST_CHECK(scalar > 1);
    BOOST_CHECK(3 > scalar);
    BOOST_CHECK(!(q1 >= q2));
    BOOST_CHECK(scalar >= 2);
    BOOST_CHECK(2 >= scalar);
    
    BOOST_CHECK_THROW(q1 < q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 <= q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 > q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 >= q3, std::runtime_error);
    
    BOOST_CHECK_THROW(q1 < scalar, std::runtime_error);
    BOOST_CHECK_THROW(scalar < q1, std::runtime_error);
    BOOST_CHECK_THROW(q1 <= scalar, std::runtime_error);
    BOOST_CHECK_THROW(scalar <= q1, std::runtime_error);
    BOOST_CHECK_THROW(q1 > scalar, std::runtime_error);
    BOOST_CHECK_THROW(scalar > q1, std::runtime_error);
    BOOST_CHECK_THROW(q1 >= scalar, std::runtime_error);
    BOOST_CHECK_THROW(scalar >= q1, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(AdditionInPlace)
{
    sycomore::Quantity q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{5, {1,0,0,0,0,0,0}};
    q1 += q2;
    BOOST_CHECK(q1 == r1);
    
    sycomore::Quantity scalar{2, {0,0,0,0,0,0,0}};
    sycomore::Quantity const r2{5, {0,0,0,0,0,0,0}};
    scalar += 3;
    BOOST_CHECK(scalar == r2);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1 += q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 += 3, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(SubtractionInPlace)
{
    sycomore::Quantity q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{-1, {1,0,0,0,0,0,0}};
    q1 -= q2;
    BOOST_CHECK(q1 == r1);
    
    sycomore::Quantity scalar{2, {0,0,0,0,0,0,0}};
    sycomore::Quantity const r2{-1, {0,0,0,0,0,0,0}};
    scalar -= 3;
    BOOST_CHECK(scalar == r2);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1 -= q3, std::runtime_error);
    BOOST_CHECK_THROW(q1 -= 3, std::runtime_error);
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

BOOST_AUTO_TEST_CASE(ConvertTo)
{
    sycomore::Quantity const q1{70, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{10, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q3{10, {0,1,0,0,0,0,0}};
    double const r = 7;
    BOOST_CHECK(q1.convert_to(q2) == r);
    BOOST_CHECK_THROW(q1.convert_to(q3), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(Double)
{
    sycomore::Quantity const scalar{3, {0,0,0,0,0,0,0}};
    sycomore::Quantity const q{3, {1,0,0,0,0,0,0}};
    BOOST_CHECK(double(scalar) == 3);
    BOOST_CHECK_THROW((double(q)), std::runtime_error);
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
    sycomore::Quantity const r1{5, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q1+q2 == r1);
    
    sycomore::Quantity const scalar{2, {0,0,0,0,0,0,0}};
    sycomore::Quantity const r2{5, {0,0,0,0,0,0,0}};
    BOOST_CHECK(scalar+3 == r2);
    BOOST_CHECK(3+scalar == r2);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1+q3, std::runtime_error);
    BOOST_CHECK_THROW(q1+3, std::runtime_error);
    BOOST_CHECK_THROW(3+q1, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(Subtraction)
{
    sycomore::Quantity const q1{2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const q2{3, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{-1, {1,0,0,0,0,0,0}};
    BOOST_CHECK(q1-q2 == r1);
    
    sycomore::Quantity const scalar{2, {0,0,0,0,0,0,0}};
    sycomore::Quantity const r2{-1, {0,0,0,0,0,0,0}};
    BOOST_CHECK(scalar-3 == r2);
    BOOST_CHECK(1-scalar == r2);

    sycomore::Quantity const q3{2, {0,1,0,0,0,0,0}};
    BOOST_CHECK_THROW(q1-q3, std::runtime_error);
    BOOST_CHECK_THROW(q1-3, std::runtime_error);
    BOOST_CHECK_THROW(1-q1, std::runtime_error);
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

BOOST_AUTO_TEST_CASE(Abs)
{
    sycomore::Quantity const q1{-9, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{9, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::abs(q1) == r1);
    
    sycomore::Quantity const q2{9, {-1,0,0,0,0,0,0}};
    sycomore::Quantity const r2{9, {-1,0,0,0,0,0,0}};
    BOOST_CHECK(std::abs(q2) == r2);
}

BOOST_AUTO_TEST_CASE(Power)
{
    sycomore::Quantity const q{9, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{3, {0.5,0,0,0,0,0,0}};
    BOOST_CHECK(std::pow(q, 0.5) == r1);
}

BOOST_AUTO_TEST_CASE(Round)
{
    sycomore::Quantity const q1{9.2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{9, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::round(q1) == r1);
    
    sycomore::Quantity const q2{-9.7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r2{-10, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::round(q2) == r2);
    
    sycomore::Quantity const q3{9, {-1.5,0,0,0,0,0,0}};
    sycomore::Quantity const r3{9, {-1.5,0,0,0,0,0,0}};
    BOOST_CHECK(std::round(q3) == r3);
}

BOOST_AUTO_TEST_CASE(Trunc)
{
    sycomore::Quantity const q1{9.2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{9, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::trunc(q1) == r1);
    
    sycomore::Quantity const q2{-9.7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r2{-9, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::trunc(q2) == r2);
    
    sycomore::Quantity const q3{9, {-1.5,0,0,0,0,0,0}};
    sycomore::Quantity const r3{9, {-1.5,0,0,0,0,0,0}};
    BOOST_CHECK(std::trunc(q3) == r3);
}

BOOST_AUTO_TEST_CASE(Floor)
{
    sycomore::Quantity const q1{9.2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{9, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::floor(q1) == r1);
    
    sycomore::Quantity const q2{-9.7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r2{-10, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::floor(q2) == r2);
    
    sycomore::Quantity const q3{9, {-1.5,0,0,0,0,0,0}};
    sycomore::Quantity const r3{9, {-1.5,0,0,0,0,0,0}};
    BOOST_CHECK(std::floor(q3) == r3);
}

BOOST_AUTO_TEST_CASE(Ceil)
{
    sycomore::Quantity const q1{9.2, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r1{10, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::ceil(q1) == r1);
    
    sycomore::Quantity const q2{-9.7, {1,0,0,0,0,0,0}};
    sycomore::Quantity const r2{-9, {1,0,0,0,0,0,0}};
    BOOST_CHECK(std::ceil(q2) == r2);
    
    sycomore::Quantity const q3{9, {-1.5,0,0,0,0,0,0}};
    sycomore::Quantity const r3{9, {-1.5,0,0,0,0,0,0}};
    BOOST_CHECK(std::ceil(q3) == r3);
}
