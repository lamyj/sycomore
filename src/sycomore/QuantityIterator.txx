#ifndef _f30849d6_016a_4c6f_9783_316ea5dfe9d0
#define _f30849d6_016a_4c6f_9783_316ea5dfe9d0

#include "QuantityIterator.h"
namespace sycomore
{

template<typename T>
QuantityIterator<T>
::QuantityIterator(T & q, bool at_end)
: _iterator(at_end ? q.magnitude.end(): q.magnitude.begin()),
    _dimensions(q.dimensions)
{
    // Nothing else
}

template<typename T>
bool
QuantityIterator<T>
::operator==(QuantityIterator<T> const & other) const
{
    return this->_iterator == other._iterator
        && this->_dimensions == other._dimensions;
}

template<typename T>
bool
QuantityIterator<T>
::operator!=(QuantityIterator<T> const & other) const
{
    return !this->operator==(other);
}

template<typename T>
QuantityConstReference
QuantityIterator<T>
::operator*() const
{
    return {*this->_iterator, this->_dimensions};
}

template<typename T>
QuantityReference
QuantityIterator<T>
::operator*()
{
    return {*this->_iterator, this->_dimensions};
}

template<typename T>
QuantityIterator<T> &
QuantityIterator<T>
::operator++()
{
    ++this->_iterator;
    return *this;
}

template<typename T>
QuantityIterator<T>
QuantityIterator<T>
::operator++(int)
{
    auto old = *this;
    this->operator++();
    return old;
}

template<typename T>
QuantityIterator<T> &
QuantityIterator<T>
::operator--()
{
    --this->_iterator;
    return *this;
}

template<typename T>
QuantityIterator<T>
QuantityIterator<T>
::operator--(int)
{
    auto const old = *this;
    this->operator--();
    return old;
}

}

#endif // _f30849d6_016a_4c6f_9783_316ea5dfe9d0
