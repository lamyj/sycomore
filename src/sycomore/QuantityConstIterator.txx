#ifndef _c45530bf_43b8_4b92_b39a_57bae62fe45d
#define _c45530bf_43b8_4b92_b39a_57bae62fe45d

#include "QuantityConstIterator.h"
namespace sycomore
{

template<typename T>
QuantityConstIterator<T>
::QuantityConstIterator(T const & q, bool at_end)
: _iterator(at_end ? q.magnitude.end(): q.magnitude.begin()),
    _dimensions(q.dimensions)
{
    // Nothing else
}

template<typename T>
bool
QuantityConstIterator<T>
::operator==(QuantityConstIterator<T> const & other) const
{
    return this->_iterator == other._iterator
        && this->_dimensions == other._dimensions;
}

template<typename T>
bool
QuantityConstIterator<T>
::operator!=(QuantityConstIterator<T> const & other) const
{
    return !this->operator==(other);
}

template<typename T>
QuantityConstReference
QuantityConstIterator<T>
::operator*() const
{
    return {*this->_iterator, this->_dimensions};
}

template<typename T>
QuantityConstIterator<T> &
QuantityConstIterator<T>
::operator++()
{
    ++this->_iterator;
    return *this;
}

template<typename T>
QuantityConstIterator<T>
QuantityConstIterator<T>
::operator++(int)
{
    auto old = *this;
    this->operator++();
    return old;
}

template<typename T>
QuantityConstIterator<T> &
QuantityConstIterator<T>
::operator--()
{
    --this->_iterator;
    return *this;
}

template<typename T>
QuantityConstIterator<T>
QuantityConstIterator<T>
::operator--(int)
{
    auto const old = *this;
    this->operator--();
    return old;
}

}

#endif // _c45530bf_43b8_4b92_b39a_57bae62fe45d
