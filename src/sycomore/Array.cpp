#include "Array.h"

#include <vector>

namespace sycomore
{

IndexGenerator
::IndexGenerator(Shape const & shape)
: _shape(shape)
{
    // Nothing else
}

IndexGenerator::const_iterator
IndexGenerator
::begin() const
{
    return const_iterator(this->_shape, false);
}

IndexGenerator::const_iterator
IndexGenerator
::end() const
{
    return const_iterator(this->_shape, true);
}

IndexGenerator::const_iterator
::const_iterator(Shape const & shape, bool is_end)
: _shape(shape), _index(this->_shape.size(), 0), _end(this->_index)
{
    std::transform(
        shape.begin(), shape.end(), this->_end.begin(),
        [](int x) { return x-1; });
    this->_next(this->_end);

    if(is_end)
    {
        this->_index = this->_end;
    }
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
    this->_next(this->_index);
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

void
IndexGenerator::const_iterator
::_next(Index & index)
{
    unsigned int i=0;
    bool increment = true;
    while(increment)
    {
        ++index[i];
        if(index[i] == this->_shape[i] && i < index.size()-1)
        {
            increment = true;
            index[i] = 0;
        }
        else
        {
            increment = false;
        }
        ++i;
    }
}

}
