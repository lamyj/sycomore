#include "GridScanner.h"

#include <algorithm>
#include <utility>

#include "sycomore/sycomore.h"

namespace sycomore
{

GridScanner
::GridScanner(Index const & origin, Shape const & shape)
: GridScanner(origin, shape, origin, shape)
{
    // Nothing else.
}

GridScanner
::GridScanner(
    Index const & origin, Shape const & shape,
    Index const & region_origin, Shape const & region_shape)
: _region_origin(region_origin), _region_end(region_origin.size())
{
    Stride stride;
    if(!shape.empty())
    {
        stride = Stride(shape.size()+1);
        stride[0] = 1;
        for(unsigned int i=0; i< shape.size(); ++i)
        {
            stride[i+1] = shape[i]*stride[i];
        }
    }

    this->_iterator.first = region_origin;
    this->_iterator.second = dot(region_origin-origin, stride);
    std::transform(
        region_origin.begin(), region_origin.end(), region_shape.begin(),
        this->_region_end.begin(), std::plus<int>());

    this->_offset_increment = Stride(stride.size());
    for(std::size_t i=0; i<this->_offset_increment.size(); ++i)
    {
        this->_offset_increment[i] = (
                (origin[i]+shape[i])-(region_origin[i]+region_shape[i])
                + region_origin[i]-origin[i]
            )*stride[i];
    }

    this->_end_index = region_origin;
    this->_end_index[this->_end_index.size()-1] =
        region_origin[this->_end_index.size()-1]
        +region_shape[this->_end_index.size()-1];
    this->_end_offset = dot(this->_end_index-origin, stride);
}

GridScanner
GridScanner
::begin() const
{
    return *this;
}

GridScanner
GridScanner
::end() const
{
    auto result(*this);
    result._iterator = std::make_pair(result._end_index, result._end_offset);
    return result;
}

GridScanner::value_type const &
GridScanner
::operator*() const
{
    return this->_iterator;
}

GridScanner::value_type const *
GridScanner
::operator->() const
{
    return &this->_iterator;
}

bool
GridScanner
::operator==(GridScanner const & other) const
{
    return this->_iterator.second == other._iterator.second;
}

bool
GridScanner
::operator!=(GridScanner const & other) const
{
    return !this->operator==(other);
}

GridScanner const &
GridScanner
::operator++()
{
    ++this->_iterator.second;

    for(std::size_t i=0, end(this->_iterator.first.size()); i<end; ++i)
    {
        ++this->_iterator.first[i];
        if(this->_iterator.first[i] == this->_region_end[i] && i!=end-1)
        {
            this->_iterator.first[i] = this->_region_origin[i];
            this->_iterator.second += this->_offset_increment[i];
        }
        else
        {
            break;
        }
    }
    return *this;
}

Index const &
GridScanner
::region_end() const
{
    return this->_region_end;
}

}
