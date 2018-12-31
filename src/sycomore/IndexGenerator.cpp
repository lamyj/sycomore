#include "IndexGenerator.h"

#include <algorithm>
#include <functional>

#include "sycomore/sycomore.h"

namespace sycomore
{

IndexGenerator
::IndexGenerator(Shape const & shape)
: IndexGenerator(Index(shape.size(), 0), shape)
{
    // Nothing else
}

IndexGenerator
::IndexGenerator(Index const & origin, Shape const & shape)
: _origin(origin), _shape(shape), _last(_origin.size()), _end(_origin.size())
{
    std::transform(
        this->_origin.begin(), this->_origin.end(), this->_shape.begin(),
        this->_last.begin(), std::plus<int>());
    this->_end = this->_origin;
    this->_end[this->_end.size()-1] = this->_last[this->_end.size()-1];
}

Index const &
IndexGenerator
::origin() const
{
    return this->_origin;
}

Shape const &
IndexGenerator
::shape() const
{
    return this->_shape;
}

Index const &
IndexGenerator
::last() const
{
    return this->_last;
}

IndexGenerator::const_iterator
IndexGenerator
::begin() const
{
    return {this->_origin, this->_shape, this->_last, this->_end, false};
}

IndexGenerator::const_iterator
IndexGenerator
::end() const
{
    return {this->_origin, this->_shape, this->_last, this->_end, true};
}

IndexGenerator::const_iterator
::const_iterator(
    Index const & origin, Shape const & shape,
    Index const & last, Index const & end, bool is_end)
: _origin(origin), _shape(shape), _last(last), _end(end),
    _index(is_end ? end : origin)
{
    // Nothing else.
}

bool
IndexGenerator::const_iterator
::operator==(const_iterator const & other) const
{
    return this->_index == other._index;
}

bool
IndexGenerator::const_iterator
::operator!=(const_iterator const & other) const
{
    return this->_index != other._index;
}

IndexGenerator::const_iterator &
IndexGenerator::const_iterator
::operator++()
{
    unsigned int i=0;
    bool increment = true;
    while(increment)
    {
        ++this->_index[i];
        if(this->_index[i] == this->_last[i] && i < this->_index.size()-1)
        {
            increment = true;
            this->_index[i] = this->_origin[i];
        }
        else
        {
            increment = false;
        }
        ++i;
    }

    return *this;
}

IndexGenerator::const_iterator &
IndexGenerator::const_iterator
::operator++(int)
{
    auto && old = *this;
    ++(*this);
    return old;
}

IndexGenerator::const_iterator::value_type const &
IndexGenerator::const_iterator
::operator*()
{
    return this->_index;
}

}
