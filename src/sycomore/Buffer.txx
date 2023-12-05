#ifndef _5de2171e_1b04_496b_aa89_b0da7fa160f2
#define _5de2171e_1b04_496b_aa89_b0da7fa160f2

#include "Buffer.h"

#include <cstddef>
#include <type_traits>
#include <xsimd/xsimd.hpp>

namespace sycomore
{

template<typename T, typename Enable>
Buffer<T, Enable>
::Buffer(size_type size)
// NOTE: new T[size] calls the constructor, but allocator.allocate(size))
// does not
: _allocator(), _data(nullptr), _capacity(0), _size(0)
{
    if(size != 0)
    {
        this->_allocate(size);
        this->_size = size;
    }
}

template<typename T, typename Enable>
Buffer<T, Enable>
::Buffer(Self const & other)
: Buffer(other._capacity)
{
    this->_size = other._size;
    for(size_type i=0; i!=this->_size; ++i)
    {
        this->_data[i] = other._data[i];
    }
}

template<typename T, typename Enable>
Buffer<T, Enable>
::Buffer(Self && other)
: _allocator(std::move(other._allocator)), _data(std::move(other._data)),
    _capacity(std::move(other._capacity)), _size(std::move(other._size))
{
    other._data = nullptr;
    other._capacity = 0;
    other._size = 0;
}

template<typename T, typename Enable>
Buffer<T, Enable>
::Buffer(std::initializer_list<T> l)
: Buffer(l.size())
{
    std::copy(l.begin(), l.end(), this->_data);
}

template<typename T, typename Enable>
Buffer<T> &
Buffer<T, Enable>
::operator=(Self const & other)
{
    if(other._size > this->_capacity)
    {
        this->_allocate(other._capacity);
    }
    
    this->_size = other._size;
    for(size_type i=0; i!=this->_size; ++i)
    {
        this->_data[i] = other._data[i];
    }
    
    return *this;
}

template<typename T, typename Enable>
Buffer<T> &
Buffer<T, Enable>
::operator=(Self && other)
{
    if(this->_data != nullptr)
    {
        this->_allocator.deallocate(this->_data, this->_size);
    }
    
    this->_allocator = std::move(other._allocator);
    this->_data = std::move(other._data);
    this->_capacity = std::move(other._capacity);
    this->_size = std::move(other._size);
    
    other._data = nullptr;
    other._capacity = 0;
    other._size = 0;
    
    return *this;
}

template<typename T, typename Enable>
Buffer<T, Enable>
::~Buffer()
{
    this->_deallocate();
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::size_type
Buffer<T, Enable>
::size() const
{
    return this->_size;
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::size_type
Buffer<T, Enable>
::capacity() const
{
    return this->_capacity;
}

template<typename T, typename Enable>
T *
Buffer<T, Enable>
::data()
{
    return this->_data;
}

template<typename T, typename Enable>
T const *
Buffer<T, Enable>
::data() const
{
    return this->_data;
}

template<typename T, typename Enable>
T const &
Buffer<T, Enable>
::operator[](std::size_t index) const
{
    return this->_data[index];
}

template<typename T, typename Enable>
T &
Buffer<T, Enable>
::operator[](std::size_t index)
{
    return this->_data[index];
}

template<typename T, typename Enable>
void
Buffer<T, Enable>
::resize(std::size_t size)
{
    if(size > this->_capacity)
    {
        auto const old_capacity = this->_capacity;
        auto old_data = this->_data;
        this->_data = nullptr;
        this->_allocate(size);
        
        std::copy(old_data, old_data+this->_size, this->_data);
        
        this->_allocator.deallocate(old_data, old_capacity);
    }
    this->_size = size;
}

template<typename T, typename Enable>
void
Buffer<T, Enable>
::resize(std::size_t size, T && x)
{
    auto const old_size = this->_size;
    this->resize(size);
    for(size_type i=old_size; i!=this->_size; ++i)
    {
        this->_data[i] = x;
    }
}

template<typename T, typename Enable>
void
Buffer<T, Enable>
::swap(Self & other)
{
    std::swap(this->_allocator, other._allocator);
    std::swap(this->_data, other._data);
    std::swap(this->_size, other._size);
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::const_iterator
Buffer<T, Enable>
::begin() const
{
    return this->_data;
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::const_iterator
Buffer<T, Enable>
::cbegin() const
{
    return this->_data;
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::const_iterator
Buffer<T, Enable>
::end() const
{
    return this->_data+this->_size;
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::const_iterator
Buffer<T, Enable>
::cend() const
{
    return this->_data+this->_size;
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::iterator
Buffer<T, Enable>
::begin()
{
    return this->_data;
}

template<typename T, typename Enable>
typename Buffer<T, Enable>::iterator
Buffer<T, Enable>
::end()
{
    return this->_data+this->_size;
}

template<typename T, typename Enable>
void
Buffer<T, Enable>
::_deallocate()
{
    if(this->_data != nullptr)
    {
        this->_allocator.deallocate(this->_data, this->_capacity);
    }
}

template<typename T, typename Enable>
void
Buffer<T, Enable>
::_allocate(std::size_t n)
{
    this->_deallocate();
    this->_data = this->_allocator.allocate(n);
    this->_capacity = n;
}

template<typename T>
void swap(Buffer<T> & b1, Buffer<T> & b2)
{
    b1.swap(b2);
}

}

#endif // _5de2171e_1b04_496b_aa89_b0da7fa160f2
