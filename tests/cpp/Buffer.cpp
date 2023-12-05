#define BOOST_TEST_MODULE Buffer
#include <boost/test/unit_test.hpp>

#include <memory>
#include <numeric>

#include "sycomore/Buffer.h"
#include "sycomore/sycomore.h"

BOOST_AUTO_TEST_CASE(EmptyPOD)
{
    sycomore::Buffer<int> b;
    BOOST_CHECK(b.size() == 0);
    BOOST_CHECK(b.capacity() == 0);
    BOOST_CHECK(b.begin() == b.end());
}

BOOST_AUTO_TEST_CASE(EmptyNonPOD)
{
    sycomore::Buffer<sycomore::Complex> b;
    BOOST_CHECK(b.size() == 0);
    BOOST_CHECK(b.capacity() == 0);
    BOOST_CHECK(b.begin() == b.end());
}

BOOST_AUTO_TEST_CASE(NonEmpty)
{
    sycomore::Buffer<sycomore::Complex> b(100);
    BOOST_CHECK(b.size() == 100);
    BOOST_CHECK(b.capacity() == 100);
    BOOST_CHECK(b.data() == &(*b.begin()));
    BOOST_CHECK(std::distance(b.begin(), b.end()) == 100);
}

BOOST_AUTO_TEST_CASE(Access)
{
    sycomore::Buffer<int> b(100);
    std::iota(b.begin(), b.end(), 0);
    
    auto it = b.begin();
    for(std::size_t i=0; i!=b.size(); ++i)
    {
        BOOST_CHECK(*it == i);
        BOOST_CHECK(b[i] == i);
        ++it;
    }
}

BOOST_AUTO_TEST_CASE(CopyConstructor)
{
    sycomore::Buffer<int> b1(100);
    std::iota(b1.begin(), b1.end(), 0);
    
    sycomore::Buffer<int> b2(b1);
    
    BOOST_CHECK(b2.size() == 100);
    BOOST_CHECK(b2.capacity() == 100);
    BOOST_CHECK(b2.data() == &(*b2.begin()));
    BOOST_CHECK(std::distance(b2.begin(), b2.end()) == 100);
    for(std::size_t i=0; i!=b2.size(); ++i)
    {
        BOOST_CHECK(b2[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(MoveConstructor)
{
    sycomore::Buffer<int> b1(100);
    std::iota(b1.begin(), b1.end(), 0);
    
    sycomore::Buffer<int> b2(std::move(b1));
    
    BOOST_CHECK(b1.size() == 0);
    BOOST_CHECK(b1.capacity() == 0);
    BOOST_CHECK(b1.data() == nullptr);
    
    BOOST_CHECK(b2.size() == 100);
    BOOST_CHECK(b2.capacity() == 100);
    BOOST_CHECK(b2.data() == &(*b2.begin()));
    BOOST_CHECK(std::distance(b2.begin(), b2.end()) == 100);
    for(std::size_t i=0; i!=b2.size(); ++i)
    {
        BOOST_CHECK(b2[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(InitConstructor)
{
    sycomore::Buffer<int> b({1,2,3});
    
    BOOST_CHECK(b.size() == 3);
    BOOST_CHECK(b.capacity() == 3);
    for(std::size_t i=0; i!=b.size(); ++i)
    {
        BOOST_CHECK(b[i] == 1+i);
    }
}

BOOST_AUTO_TEST_CASE(CopyAssignment)
{
    sycomore::Buffer<int> b1(100);
    std::iota(b1.begin(), b1.end(), 0);
    
    sycomore::Buffer<int> b2;
    b2 = b1;
    
    BOOST_CHECK(b2.size() == 100);
    BOOST_CHECK(b2.capacity() == 100);
    BOOST_CHECK(b2.data() == &(*b2.begin()));
    BOOST_CHECK(std::distance(b2.begin(), b2.end()) == 100);
    for(std::size_t i=0; i!=b2.size(); ++i)
    {
        BOOST_CHECK(b2[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(MoveAssignment)
{
    sycomore::Buffer<int> b1(100);
    std::iota(b1.begin(), b1.end(), 0);
    
    sycomore::Buffer<int> b2;
    b2 = std::move(b1);
    
    BOOST_CHECK(b1.size() == 0);
    BOOST_CHECK(b1.capacity() == 0);
    BOOST_CHECK(b1.data() == nullptr);
    
    BOOST_CHECK(b2.size() == 100);
    BOOST_CHECK(b2.capacity() == 100);
    BOOST_CHECK(b2.data() == &(*b2.begin()));
    BOOST_CHECK(std::distance(b2.begin(), b2.end()) == 100);
    
    for(std::size_t i=0; i!=b2.size(); ++i)
    {
        BOOST_CHECK(b2[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(ResizeSmaller)
{
    sycomore::Buffer<int> b(100);
    std::iota(b.begin(), b.end(), 0);
    auto const data = b.data();
    
    b.resize(50);
    
    BOOST_CHECK(b.size() == 50);
    BOOST_CHECK(b.capacity() == 100);
    BOOST_CHECK(b.data() == data);
    for(std::size_t i=0; i!=b.size(); ++i)
    {
        BOOST_CHECK(b[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(ResizeLarger)
{
    sycomore::Buffer<int> b(100);
    std::iota(b.begin(), b.end(), 0);
    auto const old_size = b.size();
    auto const old_data = b.data();
    
    b.resize(200);
    
    BOOST_CHECK(b.size() == 200);
    BOOST_CHECK(b.capacity() == 200);
    BOOST_CHECK(b.data() != old_data);
    for(std::size_t i=0; i!=old_size; ++i)
    {
        BOOST_CHECK(b[i] == i);
    }
}

BOOST_AUTO_TEST_CASE(ResizeInitLarger)
{
    sycomore::Buffer<int> b(100);
    std::iota(b.begin(), b.end(), 0);
    auto const old_size = b.size();
    auto const old_data = b.data();
    
    b.resize(200, 42);
    
    BOOST_CHECK(b.size() == 200);
    BOOST_CHECK(b.capacity() == 200);
    BOOST_CHECK(b.data() != old_data);
    for(std::size_t i=0; i!=old_size; ++i)
    {
        BOOST_CHECK(b[i] == i);
    }
    for(std::size_t i=old_size; i!=b.size(); ++i)
    {
        BOOST_CHECK(b[i] == 42);
    }
}

BOOST_AUTO_TEST_CASE(Swap)
{
    sycomore::Buffer<int> b1(100);
    std::iota(b1.begin(), b1.end(), 0);
    
    sycomore::Buffer<int> b2(10);
    std::iota(b2.begin(), b2.end(), 100);
    
    auto const b1_data = b1.data();
    auto const b2_data = b2.data();
    
    std::swap(b1, b2);
    
    BOOST_CHECK(b1.size() == 10);
    BOOST_CHECK(b1.capacity() == 10);
    BOOST_CHECK(b1.data() == b2_data);
    
    BOOST_CHECK(b2.size() == 100);
    BOOST_CHECK(b2.capacity() == 100);
    BOOST_CHECK(b2.data() == b1_data);
}
