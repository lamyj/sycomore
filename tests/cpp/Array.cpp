#define BOOST_TEST_MODULE Array
#include <boost/test/unit_test.hpp>

#include "sycomore/Array.h"

BOOST_AUTO_TEST_CASE(EmptyConstructor)
{
    sycomore::Array<int> const array;
    BOOST_TEST(array.size() == 0);
    BOOST_TEST(array.empty());
    BOOST_TEST(array.begin() == array.end());
    BOOST_TEST(array.cbegin() == array.cend());
}

BOOST_AUTO_TEST_CASE(UninitializedConstructor)
{
    sycomore::Array<int> const array(3);
    BOOST_TEST(array.size() == 3);
    BOOST_TEST(!array.empty());
    BOOST_TEST(array.begin() != array.end());
    BOOST_TEST(array.cbegin() != array.cend());
}

BOOST_AUTO_TEST_CASE(InitializedConstructor)
{
    sycomore::Array<int> const array(3, 42);
    BOOST_TEST(array.size() == 3);
    BOOST_TEST(!array.empty());
    for(size_t i=0; i<array.size(); ++i)
    {
        BOOST_TEST(array[i] == 42);
    }
}

BOOST_AUTO_TEST_CASE(CopyConstructor)
{
    sycomore::Array<size_t> array_1(3);
    for(size_t i=0; i<array_1.size(); ++i)
    {
        array_1[i] = i;
    }

    sycomore::Array<size_t> const array_2(array_1);
    BOOST_TEST(array_2.size() == array_1.size());
    BOOST_TEST(array_2.empty() == array_1.empty());
    for(size_t i=0; i<array_1.size(); ++i)
    {
        BOOST_TEST(array_2[i] == array_1[i]);
    }
}

BOOST_AUTO_TEST_CASE(InitializerListConstructor)
{
    sycomore::Array<size_t> const array{1,2,3};
    BOOST_TEST(array.size() == 3);
    BOOST_TEST(!array.empty());
    for(size_t i=0; i<array.size(); ++i)
    {
        BOOST_TEST(array[i] == i+1);
    }
}

BOOST_AUTO_TEST_CASE(MoveConstructor)
{
    sycomore::Array<size_t> array_1(3);
    for(size_t i=0; i<array_1.size(); ++i)
    {
        array_1[i] = i;
    }

    sycomore::Array<size_t> const array_2(std::move(array_1));
    BOOST_TEST(array_2.size() == 3);
    BOOST_TEST(!array_2.empty());
    for(size_t i=0; i<array_2.size(); ++i)
    {
        BOOST_TEST(array_2[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(Assignment)
{
    sycomore::Array<int> const array_1{1,2,3};
    sycomore::Array<int> array_2;
    array_2 = array_1;

    BOOST_TEST(array_2.size() == array_1.size());
    BOOST_TEST(array_2.empty() == array_1.empty());
    for(size_t i=0; i<array_1.size(); ++i)
    {
        BOOST_TEST(array_2[i] == array_1[i]);
    }
}

BOOST_AUTO_TEST_CASE(MoveAssignment)
{
    sycomore::Array<int> array_1{0,1,2};
    sycomore::Array<int> array_2;
    array_2 = array_1;

    BOOST_TEST(array_2.size() == 3);
    BOOST_TEST(!array_2.empty());
    for(size_t i=0; i<array_2.size(); ++i)
    {
        BOOST_TEST(array_2[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(Cast)
{
    sycomore::Array<unsigned int> const array{1,2,3};
    sycomore::Array<int> const e = array.astype<int>();
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 1); BOOST_TEST(e[1] == 2); BOOST_TEST(e[2] == 3);
}

BOOST_AUTO_TEST_CASE(MutatingPlusScalar)
{
    sycomore::Array<int> array{1,2,3};
    array += 4;

    BOOST_TEST(array.size() == 3);
    BOOST_TEST(array[0] == 5);
    BOOST_TEST(array[1] == 6);
    BOOST_TEST(array[2] == 7);
}

BOOST_AUTO_TEST_CASE(MutatingMinusScalar)
{
    sycomore::Array<int> array{1,2,3};
    array -= 4;

    BOOST_TEST(array.size() == 3);
    BOOST_TEST(array[0] == -3);
    BOOST_TEST(array[1] == -2);
    BOOST_TEST(array[2] == -1);
}

BOOST_AUTO_TEST_CASE(MutatingMultipliesScalar)
{
    sycomore::Array<int> array{1,2,3};
    array *= 4;

    BOOST_TEST(array.size() == 3);
    BOOST_TEST(array[0] == 4);
    BOOST_TEST(array[1] == 8);
    BOOST_TEST(array[2] == 12);
}

BOOST_AUTO_TEST_CASE(MutatingDividesScalar)
{
    sycomore::Array<float> array{1,2,3};
    array /= 4;

    BOOST_TEST(array.size() == 3);
    BOOST_TEST(array[0] == 0.25);
    BOOST_TEST(array[1] == 0.5);
    BOOST_TEST(array[2] == 0.75);
}

BOOST_AUTO_TEST_CASE(MutatingPlusArray)
{
    sycomore::Array<int> array_1{1,2,3}, array_2{4,5,6}, array_3{4,5};
    array_1 += array_2;

    BOOST_TEST(array_1.size() == 3);
    BOOST_TEST(array_1[0] == 5);
    BOOST_TEST(array_1[1] == 7);
    BOOST_TEST(array_1[2] == 9);

    BOOST_CHECK_THROW(array_1 += array_3, std::runtime_error);
    BOOST_CHECK_THROW(array_3 += array_1, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(MutatingMinusArray)
{
    sycomore::Array<int> array_1{1,2,3}, array_2{6,5,4}, array_3{4,5};
    array_1 -= array_2;

    BOOST_TEST(array_1.size() == 3);
    BOOST_TEST(array_1[0] == -5);
    BOOST_TEST(array_1[1] == -3);
    BOOST_TEST(array_1[2] == -1);

    BOOST_CHECK_THROW(array_1 -= array_3, std::runtime_error);
    BOOST_CHECK_THROW(array_3 -= array_1, std::runtime_error);
}

BOOST_AUTO_TEST_CASE(Equality)
{
    sycomore::Array<int> array_1{1,2,3}, array_2{6,5,4}, array_3{4,5};
    BOOST_TEST(array_1 == array_1);
    BOOST_TEST(array_1 != array_2);
    BOOST_TEST(array_1 != array_3);
}

BOOST_AUTO_TEST_CASE(Negate)
{
    sycomore::Array<int> a{1,2,3}; auto const e = -a;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == -1); BOOST_TEST(e[1] == -2); BOOST_TEST(e[2] == -3);
}

BOOST_AUTO_TEST_CASE(PlusScalarLeft)
{
    sycomore::Array<int> a{1,2,3}; auto const e = 4+a;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 5); BOOST_TEST(e[1] == 6); BOOST_TEST(e[2] == 7);
}

BOOST_AUTO_TEST_CASE(PlusScalarRight)
{
    sycomore::Array<int> a{1,2,3}; auto const e = a+4;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 5); BOOST_TEST(e[1] == 6); BOOST_TEST(e[2] == 7);
}

BOOST_AUTO_TEST_CASE(MinusScalar)
{
    sycomore::Array<int> a{1,2,3}; auto const e = a-4;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == -3); BOOST_TEST(e[1] == -2); BOOST_TEST(e[2] == -1);
}

BOOST_AUTO_TEST_CASE(MultipliesScalarLeft)
{
    sycomore::Array<int> a{1,2,3}; auto const e = 4*a;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 4); BOOST_TEST(e[1] == 8); BOOST_TEST(e[2] == 12);
}

BOOST_AUTO_TEST_CASE(MultipliesScalarRight)
{
    sycomore::Array<int> a{1,2,3}; auto const e = a*4;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 4); BOOST_TEST(e[1] == 8); BOOST_TEST(e[2] == 12);
}

BOOST_AUTO_TEST_CASE(DividesScalar)
{
    sycomore::Array<float> a{1,2,3}; auto const e = a/4;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 0.25); BOOST_TEST(e[1] == 0.5); BOOST_TEST(e[2] == 0.75);
}

BOOST_AUTO_TEST_CASE(PlusArraySizeMatch)
{
    sycomore::Array<int> a1{1,2,3}, a2{4,5,6}; auto const e = a1 + a2;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 5); BOOST_TEST(e[1] == 7); BOOST_TEST(e[2] == 9);
}

BOOST_AUTO_TEST_CASE(PlusArraySizeMismatch)
{
    sycomore::Array<int> a1{1,2,3}, a2{4,5}; auto const e = a1 + a2;
    BOOST_TEST(e.size() == 2);
    BOOST_TEST(e[0] == 5); BOOST_TEST(e[1] == 7);
}

BOOST_AUTO_TEST_CASE(MinusArraySizeMatch)
{
    sycomore::Array<int> a1{1,2,3}, a2{6,5,4}; auto const e = a1 - a2;
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == -5); BOOST_TEST(e[1] == -3); BOOST_TEST(e[2] == -1);
}

BOOST_AUTO_TEST_CASE(MinusArraySizeMismatch)
{
    sycomore::Array<int> a1{1,2,3}, a2{5,4}; auto const e = a1 - a2;
    BOOST_TEST(e.size() == 2);
    BOOST_TEST(e[0] == -4); BOOST_TEST(e[1] == -2);
}

BOOST_AUTO_TEST_CASE(DotSizeMatch)
{
    sycomore::Array<int> a1{1,2,3}, a2{5,4,3}; auto const x = dot(a1, a2);
    BOOST_TEST(x == 22);
}

BOOST_AUTO_TEST_CASE(DotSizeMismatch)
{
    sycomore::Array<int> a1{1,2,3}, a2{5,4}; auto const x = dot(a1, a2);
    BOOST_TEST(x == 13);
}

BOOST_AUTO_TEST_CASE(Abs)
{
    sycomore::Array<int> a{-1,2,-3}; auto const e = std::abs(a);
    BOOST_TEST(e.size() == 3);
    BOOST_TEST(e[0] == 1); BOOST_TEST(e[1] == 2); BOOST_TEST(e[2] == 3);
}
